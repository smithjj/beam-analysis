#!/usr/bin/env bash
#
# compile_wsl.sh — Legacy WSL/Linux wrapper for building the Msquared MEX files.
# Prefer compile_mex.m on all platforms.
#
# Compiles:
#   +beam/private/msquared_mex.cpp  -> msquared_mex.mexa64
#       Default: OpenMP + MATLAB's MKL-backed fft2 via mexCallMATLAB.
#       (Benchmarks showed this gives ~1.5x speedup on WSL, vs ~4x on native
#        Windows with direct FFTW. WSL2 virtualization overhead negates the
#        benefit of direct FFTW for this workload.)
#
# Usage (legacy compatibility only):
#   ./compile_wsl.sh                  # build default (MKL fft2)
#   ./compile_wsl.sh fftw             # build with MATLAB's libmwfftw3.so
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
# ---------------------------------------------------------------------------

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$REPO_DIR"

# ---------------------------------------------------------------------------
# DEFAULT build: no direct FFTW; uses mexCallMATLAB("fft2") which gets
# MATLAB's MKL.  Good enough on WSL (~1.5x peak), simpler, no library deps.
# ---------------------------------------------------------------------------
read -r -d '' CXXOPTIMFLAGS <<'EOF' || true
-O3 -march=native -ffast-math -funroll-loops -fopenmp -DNDEBUG
EOF

MEX_FLAGS=(
    "CXXOPTIMFLAGS=$CXXOPTIMFLAGS"
    "LDFLAGS=\$LDFLAGS -fopenmp"
)

# ---------------------------------------------------------------------------
# OPTIONAL FFTW build (commented out by default):
# Uncomment the block below to link MATLAB's own libmwfftw3.so for direct
# FFT access.  On native Windows this gives ~4x; on WSL2 it gives only
# ~1.5x (same as MKL), so it's kept here for benchmarking only.
# ---------------------------------------------------------------------------
# read -r -d '' CXXOPTIMFLAGS <<'EOF' || true
# -O3 -march=native -ffast-math -funroll-loops -fopenmp -DNDEBUG -DUSE_FFTW
# EOF
# MATLAB_ROOT="$(dirname "$(dirname "$(readlink -f "$(which matlab)")")")"
# FFTW_LIB="$MATLAB_ROOT/bin/glnxa64/libmwfftw3.so"
# MEX_FLAGS=(
#     "CXXOPTIMFLAGS=$CXXOPTIMFLAGS"
#     "LDFLAGS=\$LDFLAGS -fopenmp"
#     "LINKLIBS=\$LINKLIBS $FFTW_LIB"
# )

build_mex() {
    local src="+beam/private/msquared_mex.cpp"
    local out="msquared_mex"
    if [[ ! -f "$src" ]]; then
        echo "ERROR: $src not found in $REPO_DIR" >&2
        return 1
    fi
    local mode="MKL fft2 via mexCallMATLAB"
    if echo "${CXXOPTIMFLAGS:-}" | grep -q USE_FFTW; then
        mode="MATLAB libmwfftw3.so direct FFT"
    fi
    echo "==> Building $out ($mode + OpenMP)"
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
    return $missing
}

usage() {
    sed -n '2,14p' "$0"
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
        msquared_mex|fftw)   build_mex ;;
        *)
            echo "ERROR: unknown target '$t'" >&2
            echo "Valid targets: msquared_mex, fftw, clean" >&2
            exit 1
            ;;
    esac
done

echo
echo "Build complete. Run the benchmark with:"
echo "  matlab -batch 'addpath(pwd); addpath(+beam\"); examples.benchmark_msquared(); exit'"
