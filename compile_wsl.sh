#!/usr/bin/env bash
#
# compile_wsl.sh — Build the Msquared MEX files on WSL/Linux.
#
# Compiles:
#   msquared_mex.cpp  -> msquared_mex.mexa64
#       Uses OpenMP + MATLAB's MKL-backed fft2 via mexCallMATLAB.
#       (Benchmarks show MKL outperforms system FFTW3 on WSL for this workload.)
#
# Usage:
#   ./compile_wsl.sh                  # build all
#   ./compile_wsl.sh msquared_mex     # build only msquared_mex
#   ./compile_wsl.sh clean            # remove built artifacts
#
# Target CPU: tuned for 12th Gen Intel Core (Alder Lake), e.g. i7-12700H.
# The -march=native flag means the binary is NOT portable to other CPUs —
# recompile if you move it to a different machine.

set -euo pipefail

# ---------------------------------------------------------------------------
# Optimal compiler settings (validated on i7-12700H, R2026a, Ubuntu 22.04)
# ---------------------------------------------------------------------------
#   -O3                    Aggressive optimization (inlining, vectorization)
#   -march=native          Enables AVX2, FMA, ADX, all Alder Lake insns.
#                           NOT portable — see note above.
#   -ffast-math            FP reassociation, fused-multiply-add, no-NaN
#                           checks. Safe for M-squared because the math
#                           expects well-behaved field envelopes.
#   -funroll-loops         Helps tight numerical kernels (centroid sums).
#   -fopenmp               OpenMP runtime. MUST also be in LDFLAGS below,
#                           otherwise you get undefined refs to omp_get_*.
#   -DNDEBUG               Strip assert() overhead.
#
#   (FFTW3 libraries are NOT linked on WSL/Linux; benchmarks showed
#    MATLAB's MKL-backed fft2 outperforms system libfftw3 for this
#    workload.  FFTW3 is still available on Windows via compile_mex.m.)
# ---------------------------------------------------------------------------

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$REPO_DIR"

# CXXOPTIMFLAGS is picked up by `mex` for compile only.
# LDFLAGS must include -fopenmp so the linker resolves omp_get_* / GOMP_*.
# LINKLIBS picks up the FFTW libraries.
read -r -d '' CXXOPTIMFLAGS <<'EOF' || true
-O3 -march=native -ffast-math -funroll-loops -fopenmp -DNDEBUG
EOF

MEX_FLAGS=(
    "CXXOPTIMFLAGS=$CXXOPTIMFLAGS"
    "LDFLAGS=\$LDFLAGS -fopenmp"
)

build_mex() {
    local src="msquared_mex.cpp"
    local out="msquared_mex"
    if [[ ! -f "$src" ]]; then
        echo "ERROR: $src not found in $REPO_DIR" >&2
        return 1
    fi
    echo "==> Building $out (OpenMP, MATLAB fft2 fallback)"
    mex "${MEX_FLAGS[@]}" -outdir "$REPO_DIR/+beam" -output "$(basename "$out")" "$src"
    echo "    -> $REPO_DIR/+beam/$out.mexa64"
}

clean_artifacts() {
    echo "==> Cleaning build artifacts"
    rm -f +beam/msquared_mex.mexa64 +beam/msquared_mex.o
    echo "    done"
}

check_deps() {
    local missing=0
    for tool in mex g++; do
        if ! command -v "$tool" >/dev/null 2>&1; then
            echo "ERROR: '$tool' not found in PATH" >&2
            missing=1
        fi
    done
    if ! ldconfig -p 2>/dev/null | grep -q libfftw3.so; then
        echo "ERROR: libfftw3 not found. Install with:" >&2
        echo "  sudo apt install libfftw3-dev libfftw3-omp3" >&2
        missing=1
    fi
    return $missing
}

usage() {
    sed -n '2,12p' "$0"
}

# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------
if [[ "${1:-}" == "clean" ]]; then
    clean_artifacts
    exit 0
fi
if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
    usage
    exit 0
fi

check_deps

targets=("$@")
if [[ ${#targets[@]} -eq 0 ]]; then
    targets=(msquared_mex)
fi

for t in "${targets[@]}"; do
    case "$t" in
        msquared_mex)   build_mex ;;
        *)
            echo "ERROR: unknown target '$t'" >&2
            echo "Valid targets: msquared_mex, clean" >&2
            exit 1
            ;;
    esac
done

echo
echo "Build complete. Run the benchmark with:"
echo "  matlab -batch 'addpath(pwd); addpath("+beam"); examples.benchmark_msquared(); exit'"