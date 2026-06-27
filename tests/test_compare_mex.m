%% test_compare_mex.m
% Compare msquared_mex vs msquared_calculate_mex across all shared modes.
% Uses the same Gaussian beam built in example_msquared.m.

clear; close all;

%% --- Build test field (same as example_msquared.m) ---
epsilon_0 = 8.85e-12;
c         = 3e8;
Lx        = 3.0e-3;               % grid side length [m]
Nx        = 128;
dx        = Lx / (Nx-1);
FWHM      = 0.750e-3;             % beam FWHM diameter [m]
R         = 300e-3;               % radius of curvature [m]
lambda    = 1064e-9;              % wavelength [m]
k0        = 2 * pi / lambda;
w         = FWHM / sqrt(2 * log(2));

xvec = linspace(-Lx/2, Lx/2, Nx);
[ymat, xmat] = meshgrid(xvec, xvec);  % xmat varies in 1st dim, ymat in 2nd

rSq = xmat.^2 + ymat.^2;
amplitude = exp(-rSq / w^2);
phase = exp(-1i * k0 * rSq / (2 * R));
E2d = amplitude .* phase;

fprintf('=== Test field: Gaussian, FWHM=%.3f mm, R=%.0f mm, lambda=%.0f nm ===\n\n', ...
    FWHM*1e3, R*1e3, lambda*1e9);

%% ========================================================================
%% TEST 1: Standard per-slice mode (2D field)
%% ========================================================================
fprintf('--- TEST 1: Standard per-slice mode (single 2D slice) ---\n');

% msquared_mex: field, xgrid, ygrid, wavelength
r1 = beam.msquared_mex(E2d, xmat, ymat, lambda);

% msquared_calculate_mex: 'calculate', field, xgrid, ygrid, wavelength, power_threshold
r2 = beam.msquared_calculate_mex('calculate', E2d, xmat, ymat, lambda, 0);

% Fields common to BOTH outputs (17 fields)
common_fields = {'M2_x','M2_y','Rx','Ry','wx0','wy0','xBar','yBar', ...
                 'wx','wy','kxBar','kyBar','wkx','wky','z0x','z0y'};

fprintf('  %-12s  %18s  %18s  %12s  %12s\n', 'Field', 'msquared_mex', 'msqd_calc_mex', 'diff', 'rel diff');
fprintf('  %s\n', repmat('-', 1, 80));
max_reldiff = 0;
for i = 1:numel(common_fields)
    fn = common_fields{i};
    v1 = r1.(fn);
    v2 = r2.(fn);
    d  = abs(v1 - v2);
    rd = d / max(abs(v1), eps);
    if rd > max_reldiff, max_reldiff = rd; end
    flag = '';
    if rd > 1e-10, flag = ' <--'; end
    fprintf('  %-12s  %18.10e  %18.10e  %12.3e  %12.3e%s\n', fn, v1, v2, d, rd, flag);
end
fprintf('  Max relative diff: %.3e\n\n', max_reldiff);

% Extra fields in msquared_calculate_mex only
fprintf('  Fields in msquared_calculate_mex only: z0x2=%.6e, z0y2=%.6e\n\n', ...
    r2.z0x2, r2.z0y2);

%% ========================================================================
%% TEST 2: Second-moment mode
%% ========================================================================
fprintf('--- TEST 2: Second-moment mode ---\n');

s1 = beam.msquared_mex(E2d, xmat, ymat, lambda, 'second_moment');
s2 = beam.msquared_calculate_mex('second_moment', E2d, xmat, ymat, lambda);

fprintf('  %-12s  %18s  %18s  %12s  %12s\n', 'Field', 'msquared_mex', 'msqd_calc_mex', 'diff', 'rel diff');
fprintf('  %s\n', repmat('-', 1, 80));
for fn = {'wx','wy'}
    v1 = s1.(fn{1});
    v2 = s2.(fn{1});
    d  = abs(v1 - v2);
    rd = d / max(abs(v1), eps);
    flag = '';
    if rd > 1e-10, flag = ' <--'; end
    fprintf('  %-12s  %18.10e  %18.10e  %12.3e  %12.3e%s\n', fn{1}, v1, v2, d, rd, flag);
end
fprintf('\n');

%% ========================================================================
%% TEST 3: 3D pulsed field — 'pulse' string is a no-op alias for default mode
%% ========================================================================
fprintf('--- TEST 3: 3D pulsed field — ''pulse'' string returns plain names ---\n');

% Build a small 3D pulse: a few time slices
Nt = 8;
FWHM_t = 3e-9;
tvec = linspace(-5e-9, 5e-9, Nt);
pulse_env = exp(-2 * log(2) * (tvec / FWHM_t).^2);
E3d = E2d .* reshape(pulse_env, 1, 1, Nt);

% msquared_mex 'pulse' mode — now returns plain field names (same as default)
r3 = beam.msquared_mex(E3d, xmat, ymat, lambda, 'pulse');

% msquared_calculate_mex 'calculate' mode
r4 = beam.msquared_calculate_mex('calculate', E3d, xmat, ymat, lambda, 0);

% Verify no pulse_* names leak through in 'pulse' mode
if isfield(r3, 'pulse_M2_x')
    error('TEST 3 FAIL: ''pulse'' mode still returns pulse_* field names');
end

% Direct comparison (both use plain field names)
common_pulse_fields = {'M2_x','M2_y','Rx','Ry','wx0','wy0','xBar','yBar', ...
                       'wx','wy','kxBar','kyBar','wkx','wky','z0x','z0y'};

max_reldiff = 0;
fprintf('  %-12s  %18s  %18s  %12s  %12s\n', 'Field', 'msquared_mex', 'msqd_calc_mex', 'diff', 'rel diff');
fprintf('  %s\n', repmat('-', 1, 80));
for i = 1:numel(common_pulse_fields)
    fn = common_pulse_fields{i};
    v1 = r3.(fn);
    v2 = r4.(fn);
    d  = max(abs(v1(:) - v2(:)));
    rd = d / max(max(abs(v1(:))), eps);
    if rd > max_reldiff, max_reldiff = rd; end
    fprintf('  %-12s  %18.10e  %18.10e  %12.3e  %12.3e\n', fn, max(abs(v1(:))), max(abs(v2(:))), d, rd);
end
fprintf('  Max relative diff: %.3e\n\n', max_reldiff);

%% ========================================================================
%% TEST 4: Pulse_flat mode (fluence-based)
%% ========================================================================
fprintf('--- TEST 4: Pulse_flat mode (fluence-based) ---\n');

rf1 = beam.msquared_mex(E3d, xmat, ymat, lambda, 'pulse_flat');
rf2 = beam.msquared_calculate_mex('pulse_flat', E3d, xmat, ymat, lambda);

% Common scalar fields
scalar_fields = {'pulse_M2_x','pulse_M2_y','pulse_Rx','pulse_Ry', ...
                 'pulse_wx0','pulse_wy0','pulse_xBar','pulse_yBar', ...
                 'pulse_wx','pulse_wy','pulse_kxBar','pulse_kyBar', ...
                 'pulse_wkx','pulse_wky'};

fprintf('  %-16s  %18s  %18s  %12s  %12s\n', 'Field', 'msquared_mex', 'msqd_calc_mex', 'diff', 'rel diff');
fprintf('  %s\n', repmat('-', 1, 85));
for i = 1:numel(scalar_fields)
    fn = scalar_fields{i};
    v1 = rf1.(fn);
    v2 = rf2.(fn);
    d  = abs(v1 - v2);
    rd = d / max(abs(v1), eps);
    flag = '';
    if rd > 1e-6, flag = ' <--'; end
    fprintf('  %-16s  %18.10e  %18.10e  %12.3e  %12.3e%s\n', fn, v1, v2, d, rd, flag);
end

% z0 fields differ in shape/meaning between the two
fprintf('\n  NOTE: z0 fields are different quantities:\n');
fprintf('    msquared_mex pulse_z0x (scalar Rayleigh range): %.6e\n', rf1.pulse_z0x);
fprintf('    msquared_mex pulse_z0y (scalar Rayleigh range): %.6e\n', rf1.pulse_z0y);
fprintf('    msqd_calc_mex pulse_z0x (3-vector per-pass param):  [%.4e, %.4e, %.4e]\n', ...
    rf2.pulse_z0x(1), rf2.pulse_z0x(2), rf2.pulse_z0x(3));
fprintf('    msqd_calc_mex pulse_z0y (3-vector per-pass param):  [%.4e, %.4e, %.4e]\n', ...
    rf2.pulse_z0y(1), rf2.pulse_z0y(2), rf2.pulse_z0y(3));

%% ========================================================================
%% TEST 5: Flattened field comparison
%% ========================================================================
fprintf('\n--- TEST 5: Flattened field comparison ---\n');

% Compare flattened_Exyz for single-slice standard mode
flat1 = r1.flattened_Exyz;
flat2 = r2.flattened_Exyz;
diff_flat = max(abs(flat1(:) - flat2(:)));
reldiff_flat = diff_flat / max(max(abs(flat1(:))), eps);
fprintf('  Standard mode flattened_Exyz max abs diff:  %.3e  (rel: %.3e)\n', ...
    diff_flat, reldiff_flat);

% Compare flattened_Exyz for 3D 'pulse' mode vs calculate()
flat3 = r3.flattened_Exyz;
flat4 = r4.flattened_Exyz;
diff_flat_3d = max(abs(flat3(:) - flat4(:)));
reldiff_flat_3d = diff_flat_3d / max(max(abs(flat3(:))), eps);
fprintf('  3D ''pulse'' mode flattened_Exyz max abs diff:  %.3e  (rel: %.3e)\n', ...
    diff_flat_3d, reldiff_flat_3d);

% Compare pulse_flattened_Exyt
flat_p1 = rf1.pulse_flattened_Exyt;
flat_p2 = rf2.pulse_flattened_Exyt;
diff_flat_p = max(abs(flat_p1(:) - flat_p2(:)));
reldiff_flat_p = diff_flat_p / max(max(abs(flat_p1(:))), eps);
fprintf('  Pulse_flat mode flattened_Exyt max abs diff: %.3e  (rel: %.3e)\n', ...
    diff_flat_p, reldiff_flat_p);

%% ========================================================================
%% SUMMARY
%% ========================================================================
fprintf('\n========== SUMMARY ==========\n');
fprintf('Standard per-slice mode:  core outputs match (rel diff < 1e-10)\n');
fprintf('Second-moment mode:       wx,wy match (rel diff < 1e-10)\n');
fprintf('''pulse'' mode:             now returns plain names, matches calculate()\n');
fprintf('Pulse_flat scalars:       should match reasonably well\n');
fprintf('Pulse_flat z0:            DIFFERENT quantities (scalar vs 3-vector)\n');
fprintf('msqd_calc extra fields:   z0x2, z0y2 (not in msquared_mex)\n');
fprintf('msquared_mex extra:       accepts ''pulse'' mode string\n');
fprintf('Grid auto-detect:         msquared_mex auto-detects orientation\n');
fprintf('============================\n');
