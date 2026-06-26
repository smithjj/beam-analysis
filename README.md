# M² Beam Quality Calculator (MEX + MATLAB Class)

This directory contains a drop-in MEX replacement for `Msquared_class.calculate()` and `Msquared_class.calculate_pulse_flat()`, plus example scripts and regression tests.

## Files

| File | Description |
|------|-------------|
| `Msquared_class.m` | Original MATLAB reference implementation (Arlee Smith / AS-Photonics). |
| `msquared_mex.cpp` | C++ MEX source implementing the same algorithms. |
| `msquared_mex.mexw64` | Compiled MEX binary (Windows x64). |
| `example_msquared.m` | CW Gaussian beam example. |
| `example_msquared_pulse.m` | Pulsed Gaussian beam example (per-slice + fluence-based). |
| `test_msquared.m` | Core regression tests (19 tests). |
| `test_msquared_extended.m` | Extended tests: displacements, curvature, super-Gaussian, noise, edge cases (27 tests). |

## Building the MEX

MATLAB must have a supported C++ compiler. If `mex -setup` only finds MinGW, but Visual Studio 2022 Community is installed, force the MSVC option file:

```matlab
mex -f 'C:\Program Files\MATLAB\R2026a\bin\win64\mexopts\msvcpp2022.xml' msquared_mex.cpp
```

> **Note:** MSVC does not support GNU `sincos()`, so `msquared_mex.cpp` uses `std::sin`/`std::cos`.

## Usage

### CW / per-slice analysis

```matlab
results = msquared_mex(field, xgrid, ygrid, wavelength);
```

- `field` — complex `Nx × Ny` array, or `Nx × Ny × Nt` array for per-slice analysis.
- `xgrid`, `ygrid` — vectors or 2-D coordinate grids (m).
- `wavelength` — wavelength (m).

Returned fields include `M2_x`, `M2_y`, `Rx`, `Ry`, `wx`, `wy`, `wx0`, `wy0`, `z0x`, `z0y`, centroids, k-space widths, and `flattened_Exyz`.

### Second-moment only

```matlab
w = msquared_mex(field, xgrid, ygrid, wavelength, 'second_moment');
```

Fast path that computes only the spatial second-moment widths. Returns a struct with fields `wx` and `wy` (`Nt × 1` arrays), without performing FFTs, curvature removal, or M² calculations.

### Pulsed per-slice analysis (no-op alias)

```matlab
results = msquared_mex(field, xgrid, ygrid, wavelength, 'pulse');
```

Identical to the default per-slice analysis above. The string `'pulse'` is accepted as a convenience alias, and output field names are **not** prefixed with `pulse_` — they match `Msquared_class.calculate()` (`M2_x`, `wx`, `z0x`, `flattened_Exyz`, etc.).

### Fluence-based pulse analysis

```matlab
results = msquared_mex(field, xgrid, ygrid, wavelength, 'pulse_flat');
% or equivalently:
results = msquared_mex(field, xgrid, ygrid, wavelength, 'calculate_pulse_flat');
```

This implements `Msquared_class.calculate_pulse_flat()`: it time-integrates the fluence first, then computes a **single scalar** M² for the whole pulse. Output field names are prefixed with `pulse_`, but they are scalars (not arrays), except for `pulse_flattened_Exyt` which remains `Nx*Ny × Nt`.

## Running the Tests

From MATLAB:

```matlab
runtests('test_msquared')
runtests('test_msquared_extended')
```

Or run both and summarize:

```matlab
clear classes;
r1 = runtests('test_msquared');
r2 = runtests('test_msquared_extended');
fprintf('Total: %d passed, %d failed, %d incomplete\n', ...
    sum([r1.Passed; r2.Passed]), ...
    sum([r1.Failed; r2.Failed]), ...
    sum([r1.Incomplete; r2.Incomplete]));
```

> **Important:** `clear classes` is needed after editing a class-definition test file, because MATLAB caches class definitions.

## Test Coverage

### Core tests (`test_msquared.m`) — 19 tests

| Test | What it checks |
|------|---------------|
| `testIdealGaussianCW_MEX` / `_class` | Ideal Gaussian CW beam → M² = 1, wx = w. |
| `testCwMEXvsClass` | MEX and MATLAB class agree to high precision for CW. |
| `testRadiusOfCurvature_MEX` | Reported \|R\| matches input radius. |
| `testCentroidAtOrigin` | Centred beam has centroid at (0, 0). |
| `testZeroField` | Zero field does not crash and returns M² = 0. |
| `testPulsePerSlice_MEX` / `_class` / `_MEXvsClass` | Per-slice pulse analysis returns M² = 1 per slice; MEX/class agree. |
| `testPulseFlat_IdealGaussian` / `_MEXvsClass` | Fluence-based analysis returns M² = 1 and matches the MATLAB class. |
| `testPulseFlat_vs_PerSliceAverage` | Uniform pulse: fluence-based result matches per-slice result. |
| `testM2_gt_1_superGaussian` | Super-Gaussian beam gives M² > 1; MEX/class agree. |
| `testZeroField_pulse_flat` | Zero 3-D field in `pulse_flat` mode does not crash. |
| `testOutputStructFields*` | All expected output fields exist with correct sizes/scalars. |
| `testVectorGridInput` | Vector vs matrix grid inputs give identical results. |
| `testWavelengthVaries` | Same M² for different wavelengths, but different Rayleigh ranges. |

### Extended tests (`test_msquared_extended.m`) — 27 tests

| Category | Tests |
|----------|-------|
| **Transverse displacements** | `testDisplacedX`, `testDisplacedY`, `testDisplacedXY`, `testDisplacedNearEdge` |
| **Radii of curvature** | `testRadius50mm`, `testRadius100mm`, `testRadius500mm`, `testRadius1000mm`, `testRadiusInf`, `testCurvatureSignConvergence` |
| **Super-Gaussian shapes** | `testSuperGaussianOrder4/6/8`, `testSuperGaussianM2IncreasesWithOrder` |
| **Noisy / messy beams** | `testAmplitudeNoise5percent`, `testAmplitudeNoise20percent`, `testPhaseNoise`, `testMessyBeamManyDefects`, `testClippedGaussian` |
| **Astigmatism & tilt** | `testAstigmaticBeam`, `testTiltedWavefrontX`, `testTiltedWavefrontY` |
| **Edge cases** | `testAllRealFieldDoesNotCrash`, `testAllRealFieldGivesWrongPhysics`, `testVerySmallFieldValues`, `testSingleTimeSliceFallsBack`, `testPulseFlatWithDisplacedGaussian` |

## Known Behaviours / Caveats

1. **All-real input fields do not crash, but are not physically valid.**
   A purely real field (e.g. `cos(k r² / 2R)` instead of `exp(-i k r² / 2R)`) will run and return numbers, but the result is not a meaningful propagating Gaussian beam. Example output for an all-real Gaussian: `M2_x ≈ 3.9`, `Rx = Inf`.

2. **Finite grid truncation.**
   Second-moment M² is sensitive to truncation. Displacing a beam near the grid edge or using a very tight radius of curvature can raise M² above 1 even for an otherwise ideal Gaussian. The tests account for this by checking agreement between MEX and class and by verifying values remain finite and reasonable, rather than demanding exactly M² = 1 in those regimes.

3. **MEX vs class field-name convention.**
   - `Msquared_class.calculate()` always returns fields named `M2_x`, `M2_y`, `wx`, etc., regardless of whether the input is 2-D or 3-D.
   - `msquared_mex` returns `M2_x`, `M2_y`, … for 2-D, unflagged 3-D, or 3-D input with `'pulse'`.
   - `'pulse_flat'` uses the same `pulse_` prefixed names as `Msquared_class.calculate_pulse_flat()`. Note that `msquared_mex` computes a scalar Rayleigh range there, while the class returns a 3-element per-pass propagation parameter.
   - `msquared_calculate_mex` mirrors the class exactly and includes extra fields `z0x2` and `z0y2` that `msquared_mex` omits.

## Current Status

```
test_msquared:          19 passed, 0 failed, 0 incomplete
test_msquared_extended: 27 passed, 0 failed, 0 incomplete
Total:                  46 passed, 0 failed, 0 incomplete
```
