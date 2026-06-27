# Msquared Toolbox — M² Beam Quality Calculator

MATLAB toolbox for calculating the M-squared beam-quality factor and related
beam quantities for 2-D (x,y) and 3-D (x,y,t) electric-field arrays.  It
includes a reference MATLAB class and compiled MEX drop-in accelerators.

## Installation

Install from the `.mltbx` file in MATLAB:

```matlab
matlab.addons.install('MsquaredToolbox.mltbx')
```

Or add the repository root to your MATLAB path:

```matlab
addpath('C:\path\to\beam-analysis')
```

The functions live in the `+beam` namespace and are accessed as
`beam.Msquared`, `beam.msquared_mex`, and `beam.msquared_calculate_mex`.

## Repository Layout

| Path | Description |
|------|-------------|
| `+beam/Msquared.m` | Reference MATLAB class implementation. |
| `+beam/msquared_mex.*` | Compiled MEX accelerators (`.mexw64`, `.mexa64`). |
| `msquared_mex.cpp` | C++ MEX source. |
| `msquared_calculate_mex.cpp` | Lower-level C++ MEX source. |
| `examples/` | Example scripts and benchmarking. |
| `tests/test_compare_mex.m` | MEX-vs-class regression tests. |
| `compile_mex.m` / `compile_wsl.sh` | MEX build scripts for Windows and WSL/Linux. |
| `ToolboxOptions.m` | Toolbox packaging configuration. |
| `package_toolbox.m` | One-click build of the installable `.mltbx`. |

## Building the MEX

### Windows

```matlab
compile_mex
```

### WSL/Linux

```bash
./compile_wsl.sh
```

Both scripts place the compiled binaries into `+beam/` so they are
automatically available through the namespace.

## Usage

### CW / per-slice analysis

```matlab
results = beam.msquared_mex(field, xgrid, ygrid, wavelength);
```

- `field` — complex `Nx × Ny` array, or `Nx × Ny × Nt` array for per-slice analysis.
- `xgrid`, `ygrid` — vectors or 2-D coordinate grids (m).
- `wavelength` — wavelength (m).

Returned fields include `M2_x`, `M2_y`, `Rx`, `Ry`, `wx`, `wy`, `wx0`, `wy0`, `z0x`, `z0y`, centroids, k-space widths, and `flattened_Exyz`.

### Using the class

```matlab
m2 = beam.Msquared();
m2.xgrid = xmat;
m2.ygrid = ymat;
m2.field = E;
results = m2.calculate();

% Fluence-based pulse analysis
results_flat = m2.calculate_pulse_flat();
```

### Second-moment only

```matlab
w = beam.msquared_mex(field, xgrid, ygrid, wavelength, 'second_moment');
```

Fast path that computes only the spatial second-moment widths. Returns a struct with fields `wx` and `wy` (`Nt × 1` arrays), without performing FFTs, curvature removal, or M² calculations.

### Pulsed per-slice analysis (no-op alias)

```matlab
results = beam.msquared_mex(field, xgrid, ygrid, wavelength, 'pulse');
```

Identical to the default per-slice analysis. The string `'pulse'` is accepted as a convenience alias; output field names are **not** prefixed with `pulse_` — they match `beam.Msquared().calculate()` (`M2_x`, `wx`, `z0x`, `flattened_Exyz`, etc.).

### Fluence-based pulse analysis

```matlab
results = beam.msquared_mex(field, xgrid, ygrid, wavelength, 'pulse_flat');
% or equivalently:
results = beam.msquared_mex(field, xgrid, ygrid, wavelength, 'calculate_pulse_flat');
```

This implements `beam.Msquared().calculate_pulse_flat()`: it time-integrates the fluence first, then computes a **single scalar** M² for the whole pulse. Output field names are prefixed with `pulse_`, but they are scalars (not arrays), except for `pulse_flattened_Exyt` which remains `Nx*Ny × Nt`.

## Running the Tests

From MATLAB:

```matlab
runtests('tests.test_compare_mex')
```

## MEX vs Class Field-Name Convention

- `beam.Msquared().calculate()` always returns fields named `M2_x`, `M2_y`, `wx`, etc., regardless of whether the input is 2-D or 3-D.
- `beam.msquared_mex` returns `M2_x`, `M2_y`, … for 2-D, unflagged 3-D, or 3-D input with `'pulse'`.
- `'pulse_flat'` uses the same `pulse_` prefixed names as `beam.Msquared().calculate_pulse_flat()`.
- `beam.msquared_calculate_mex` mirrors the class exactly and includes extra fields `z0x2` and `z0y2` that `beam.msquared_mex` omits.

## Known Behaviours / Caveats

1. **All-real input fields do not crash, but are not physically valid.**
   A purely real field will run and return numbers, but the result is not a meaningful propagating Gaussian beam.

2. **Finite grid truncation.**
   Second-moment M² is sensitive to truncation. Displacing a beam near the grid edge or using a very tight radius of curvature can raise M² above 1 even for an otherwise ideal Gaussian.

3. **Compiled binaries are platform-specific.**
   The `.mexw64` files run on Windows; the `.mexa64` files run on Linux/WSL. Re-run the appropriate compile script on a new platform.

## Package the Toolbox

```matlab
package_toolbox
```

Produces `MsquaredToolbox.mltbx`, which can be installed via MATLAB’s Add-On Manager.
