%% Example2: Pulsed M^2 measurement using msquared_mex and Msquared
%
% Same beam/pulse construction as example_msquared_pulse.m, but computes
% parameters with both the MEX function and the native MATLAB class and
% prints the two sets of numbers side-by-side for easy comparison.
%
% Grid:       3 mm x 3 mm, 128 x 128 points
% Beam:       FWHM diameter = 0.750 mm,  R = 300 mm
% Wavelength: 1064 nm
% Pulse:      FWHM = 3 ns,  t = -10 .. +10 ns, 129 pts, 1 mJ energy

%% Physical constants
epsilon_0 = 8.854187817e-12;    % vacuum permittivity [F/m]
c         = 299792458;           % speed of light [m/s]

%% Spatial parameters (SI units: meters)
Lx      = 3.0e-3;                % grid side length [m]
Nx      = 128;                   % points per side
dx      = Lx / (Nx-1);           % step size
Ly      = Lx;
Ny      = 128;
dy      = dx;

FWHM    = 0.750e-3;              % beam FWHM diameter [m]
R       = 300e-3;                % radius of curvature [m]
lambda  = 1064e-9;               % wavelength [m]

k0      = 2 * pi / lambda;       % wavenumber [rad/m]

fprintf('=== Spatial grid ===\n');
fprintf('Grid:  %d x %d points over %.2f x %.2f mm\n', Nx, Nx, Lx*1e3, Lx*1e3);
fprintf('dx:    %.3f um\n', dx*1e6);
fprintf('FWHM:  %.3f mm  |  R: %.1f mm  |  lambda: %.0f nm\n', FWHM*1e3, R*1e3, lambda*1e9);

%% Build spatial grid
%
%  msquared_mex expects columns (ix) to be the 2nd dimension and rows (iy)
%  the 1st.  The MEX convention is:
%    xgrid varies in the 1st dim (rows),   xmat(iy,ix) = xgrid(iy)
%    ygrid varies in the 2nd dim (columns), ymat(iy,ix) = ygrid(ix)
%
%  We use ndgrid to build 2-D grid matrices matching this convention.
xv = linspace(-Lx/2, Lx/2, Nx);
yv = linspace(-Ly/2, Ly/2, Ny);

% ndgrid: X(iy,ix) = xv(iy),  Y(iy,ix) = yv(ix)
[X2d, Y2d] = ndgrid(xv, xv);

%% Construct continuous-wave (CW) Gaussian field
%
%  E(x,y) = exp(-r^2/w^2) * exp(-i*k*r^2/(2R))
%
%  w : 1/e^2 intensity radius
%      FWHM = w * sqrt(2*ln2)  -->  w = FWHM / sqrt(2*ln2)

w   = FWHM / sqrt(2 * log(2));     % 1/e^2 radius [m]
rSq = X2d.^2 + Y2d.^2;             % squared radius [m^2]

E_cw = exp(-rSq / w^2) .* exp(-1i * k0 * rSq / (2 * R));

fprintf('\nBeam 1/e^2 radius w = %.4f mm\n', w*1e3);

%% Temporal parameters
t_min  = -10e-9;                   % s
t_max  =  10e-9;                   % s
Nt     = 129;                       % number of time samples
tvec   = linspace(t_min, t_max, Nt);
dt     = tvec(2) - tvec(1);        % time step [s]

fprintf('\n=== Temporal grid ===\n');
fprintf('t:     %.1f .. %.1f ns,  %d points,  dt = %.3f ns\n', ...
    t_min*1e9, t_max*1e9, Nt, dt*1e9);

FWHM_t      = 3e-9;                % pulse FWHM duration [s]
pulse_energy = 1e-3;               % total pulse energy [J]
fprintf('Pulse: FWHM = %.1f ns,  Energy = %.1f mJ\n', FWHM_t*1e9, pulse_energy*1e3);

%% Gaussian temporal pulse envelope
%
%  Intensity:  I(t)       = I0 * exp(-4*ln2 * t^2 / FWHM_t^2)
%  Field env:  A(t)/A0    = sqrt(I(t)/I0) = exp(-2*ln2 * t^2 / FWHM_t^2)
%
pulse_env = exp(-2 * log(2) * (tvec / FWHM_t).^2);
fprintf('Pulse envelope peak = %.4f (t=0)\n', max(pulse_env));

%% Build 3-D pulsed field: E_pulsed(x,y,t) = E_cw(x,y) * A(t)
%
%  Shape: Nx x Ny x Nt  (here Ny == Nx)
%
E_pulsed = E_cw .* reshape(pulse_env, 1, 1, Nt);   % implicit expansion

%% Normalise to 1 mJ total pulse energy
%
%  Energy = 0.5 * epsilon_0 * c * sum(|E|^2) * dx * dy * dt   [J]
%
energy = 0.5 * epsilon_0 * c * sum(abs(E_pulsed(:)).^2) * dx * dy * dt;
E_pulsed = E_pulsed * sqrt(pulse_energy / energy);

energy_check = 0.5 * epsilon_0 * c * sum(abs(E_pulsed(:)).^2) * dx * dy * dt;
fprintf('\nPulse energy: %.3f mJ (target %.0f mJ)\n', energy_check*1e3, pulse_energy*1e3);

%% ========================================================================
%% Per-time-slice analysis
%% ========================================================================

% --- MEX per-slice mode ---
fprintf('\nRunning msquared_mex (instantaneous mode) ...\n');
tic;
results_mex = msquared_auto(E_pulsed, xv, yv, lambda, 'pulse');
t_mex = toc;

% --- Msquared per-slice mode ---
fprintf('Running Msquared.calculate() ...\n');
m2 = beam.Msquared();
m2.field      = E_pulsed;
m2.xgrid      = X2d;          % 2-D ndgrid arrays (X varies in dim 1, Y in dim 2)
m2.ygrid      = Y2d;
m2.wavelength = lambda;
tic;
results_class = m2.calculate();
t_class = toc;

% Centre time slice (closest to t = 0)
[~, t0_idx] = min(abs(tvec));

fprintf('\n========== Per-slice beam analysis @ t = %.2f ns (MEX vs class) ==========\n', ...
    tvec(t0_idx)*1e9);
fprintf('Quantity          MEX                MATLAB class       Units\n');
fprintf('%s\n', repmat('-', 1, 70));
fprintf('M^2_x:            %.4f             %.4f\n',    results_mex.M2_x(t0_idx),    results_class.M2_x(t0_idx));
fprintf('M^2_y:            %.4f             %.4f\n',    results_mex.M2_y(t0_idx),    results_class.M2_y(t0_idx));
fprintf('wx (1/e^2):       %.4f mm          %.4f mm\n', results_mex.wx(t0_idx)*1e3,  results_class.wx(t0_idx)*1e3);
fprintf('wy (1/e^2):       %.4f mm          %.4f mm\n', results_mex.wy(t0_idx)*1e3,  results_class.wy(t0_idx)*1e3);
fprintf('wx0 (waist):      %.4f mm          %.4f mm\n', results_mex.wx0(t0_idx)*1e3, results_class.wx0(t0_idx)*1e3);
fprintf('wy0 (waist):      %.4f mm          %.4f mm\n', results_mex.wy0(t0_idx)*1e3, results_class.wy0(t0_idx)*1e3);
fprintf('Rx:               %+.1f mm         %+.1f mm\n',results_mex.Rx(t0_idx)*1e3,   results_class.Rx(t0_idx)*1e3);
fprintf('Ry:               %+.1f mm         %+.1f mm\n',results_mex.Ry(t0_idx)*1e3,   results_class.Ry(t0_idx)*1e3);
fprintf('z0x (Rayleigh):   %.2f mm          %.2f mm\n', results_mex.z0x(t0_idx)*1e3, results_class.z0x(t0_idx)*1e3);
fprintf('z0y (Rayleigh):   %.2f mm          %.2f mm\n', results_mex.z0y(t0_idx)*1e3, results_class.z0y(t0_idx)*1e3);
fprintf('xBar centroid:    %.6f mm        %.6f mm\n', results_mex.xBar(t0_idx)*1e3, results_class.xBar(t0_idx)*1e3);
fprintf('yBar centroid:    %.6f mm        %.6f mm\n', results_mex.yBar(t0_idx)*1e3, results_class.yBar(t0_idx)*1e3);
fprintf('kxBar centroid:   %.6f rad/mm     %.6f rad/mm\n', results_mex.kxBar(t0_idx)*1e-3, results_class.kxBar(t0_idx)*1e-3);
fprintf('kyBar centroid:   %.6f rad/mm     %.6f rad/mm\n', results_mex.kyBar(t0_idx)*1e-3, results_class.kyBar(t0_idx)*1e-3);
fprintf('wkx (angular):    %.6f rad/mm     %.6f rad/mm\n', results_mex.wkx(t0_idx)*1e-3, results_class.wkx(t0_idx)*1e-3);
fprintf('wky (angular):    %.6f rad/mm     %.6f rad/mm\n', results_mex.wky(t0_idx)*1e-3, results_class.wky(t0_idx)*1e-3);
fprintf('=========================================================================\n');
fprintf('MEX time:   %.4f s\n', t_mex);
fprintf('Class time: %.4f s\n', t_class);

%% ========================================================================
%% Fluence-based pulse analysis
%% ========================================================================

% --- MEX pulse_flat mode ---
fprintf('\nRunning msquared_mex (fluence-based pulse_flat mode) ...\n');
tic;
results_flat_mex = msquared_auto(E_pulsed, X2d, Y2d, lambda, 'pulse_flat');
t_flat_mex = toc;

% --- Msquared pulse_flat mode ---
fprintf('Running Msquared.calculate_pulse_flat() ...\n');
tic;
results_flat_class = m2.calculate_pulse_flat();
t_flat_class = toc;

fprintf('\n========== Fluence-based pulse analysis (MEX vs class) ==========\n');
fprintf('Quantity            MEX                MATLAB class       Units\n');
fprintf('%s\n', repmat('-', 1, 70));
fprintf('pulse_M2_x:        %.4f             %.4f\n',      results_flat_mex.pulse_M2_x,      results_flat_class.pulse_M2_x);
fprintf('pulse_M2_y:        %.4f             %.4f\n',      results_flat_mex.pulse_M2_y,      results_flat_class.pulse_M2_y);
fprintf('pulse_wx:          %.4f mm          %.4f mm\n',   results_flat_mex.pulse_wx*1e3,    results_flat_class.pulse_wx*1e3);
fprintf('pulse_wy:          %.4f mm          %.4f mm\n',   results_flat_mex.pulse_wy*1e3,    results_flat_class.pulse_wy*1e3);
fprintf('pulse_wx0:         %.4f mm          %.4f mm\n',   results_flat_mex.pulse_wx0*1e3,   results_flat_class.pulse_wx0*1e3);
fprintf('pulse_wy0:         %.4f mm          %.4f mm\n',   results_flat_mex.pulse_wy0*1e3,   results_flat_class.pulse_wy0*1e3);
fprintf('pulse_Rx:          %+.1f mm         %+.1f mm\n',  results_flat_mex.pulse_Rx*1e3,    results_flat_class.pulse_Rx*1e3);
fprintf('pulse_Ry:          %+.1f mm         %+.1f mm\n',  results_flat_mex.pulse_Ry*1e3,    results_flat_class.pulse_Ry*1e3);
fprintf('pulse_xBar:        %.6f mm        %.6f mm\n',   results_flat_mex.pulse_xBar*1e3,  results_flat_class.pulse_xBar*1e3);
fprintf('pulse_yBar:        %.6f mm        %.6f mm\n',   results_flat_mex.pulse_yBar*1e3,  results_flat_class.pulse_yBar*1e3);
fprintf('pulse_kxBar:       %.6f rad/mm     %.6f rad/mm\n', results_flat_mex.pulse_kxBar*1e-3, results_flat_class.pulse_kxBar*1e-3);
fprintf('pulse_kyBar:       %.6f rad/mm     %.6f rad/mm\n', results_flat_mex.pulse_kyBar*1e-3, results_flat_class.pulse_kyBar*1e-3);
fprintf('pulse_wkx:         %.6f rad/mm     %.6f rad/mm\n', results_flat_mex.pulse_wkx*1e-3,  results_flat_class.pulse_wkx*1e-3);
fprintf('pulse_wky:         %.6f rad/mm     %.6f rad/mm\n', results_flat_mex.pulse_wky*1e-3,  results_flat_class.pulse_wky*1e-3);
% z0 fields differ in shape/meaning between the two implementations
fprintf('pulse_z0x:         %.4f mm          %.4f mm*\n', ...
    results_flat_mex.pulse_z0x*1e3, results_flat_class.pulse_z0x(1)*1e3);
fprintf('pulse_z0y:         %.4f mm          %.4f mm*\n', ...
    results_flat_mex.pulse_z0y*1e3, results_flat_class.pulse_z0y(1)*1e3);
fprintf('=================================================================\n');
fprintf('MEX time:   %.4f s\n', t_flat_mex);
fprintf('Class time: %.4f s\n', t_flat_class);

%% Visualization

figure('Position', [50 50 1400 800]);

% --- 1. Fluence map (time-integrated intensity) ---
fluency = 0.5 * epsilon_0 * c * trapz(abs(E_pulsed).^2, 3) * dt;  % J/m^2

subplot(2, 3, 1);
imagesc(xv*1e6, xv*1e6, fluency);
axis image; axis xy;
xlabel('x (\mum)'); ylabel('y (\mum)');
title('Fluence [J/m^2]');
colorbar; colormap(gca, 'hot');

% --- 2. Temporal profile at beam centre ---
[~, cx] = min(abs(xv));
I_center = 0.5 * epsilon_0 * c * squeeze(abs(E_pulsed(cx, cx, :))).^2;  % W/m^2
I_display = I_center * 1e-4;   % convert to W/cm^2 for plotting

% --- FWHM measurement with interpolation ---
I_peak   = max(I_display);
half_max = I_peak / 2;
above_half = I_display >= half_max;
idx_left  = find(above_half, 1, 'first');
idx_right = find(above_half, 1, 'last');

% Linear interpolation for sub-sample accurate FWHM
if idx_left > 1
    t_left = interp1(I_display(idx_left-1:idx_left), ...
                     tvec(idx_left-1:idx_left), half_max);
else
    t_left = tvec(idx_left);
end
if idx_right < length(tvec)
    t_right = interp1(I_display(idx_right:idx_right+1), ...
                      tvec(idx_right:idx_right+1), half_max);
else
    t_right = tvec(idx_right);
end
t_fwhm = (t_right - t_left) * 1e9;

subplot(2, 3, 2);
plot(tvec*1e9, I_display, 'b-', 'LineWidth', 1.5); hold on;
yline(half_max, 'r--', 'LineWidth', 1);
xline([t_left t_right]*1e9, 'g--', 'LineWidth', 1);
xlabel('Time (ns)'); ylabel('Intensity (W/cm^2)');
title(sprintf('Temp. profile FWHM = %.2f ns', t_fwhm));
legend('Intensity', 'Half-max', 'FWHM', 'Location', 'north');
grid on;

% --- 3. M^2 vs time (should be flat for an ideal Gaussian) ---
subplot(2, 3, 3);
plot(tvec*1e9, results_mex.M2_x, 'b.-', ...
    'LineWidth', 1, 'MarkerSize', 8); hold on;
plot(tvec*1e9, results_class.M2_x, 'r.--', ...
    'LineWidth', 1, 'MarkerSize', 8);
plot(tvec*1e9, results_mex.M2_y, 'c.-', ...
    'LineWidth', 1, 'MarkerSize', 8);
plot(tvec*1e9, results_class.M2_y, 'm.--', ...
    'LineWidth', 1, 'MarkerSize', 8);
xlabel('Time (ns)'); ylabel('M^2');
title('M^2 vs time');
legend('MEX M^2_x', 'Class M^2_x', 'MEX M^2_y', 'Class M^2_y', 'Location', 'best');
ylim([min(results_mex.M2_x)*0.999, max(results_mex.M2_x)*1.001]);
grid on;

% --- 4. Beam width vs time ---
subplot(2, 3, 4);
plot(tvec*1e9, results_mex.wx*1e6, 'b.-', ...
    'LineWidth', 1, 'MarkerSize', 8); hold on;
plot(tvec*1e9, results_class.wx*1e6, 'r.--', ...
    'LineWidth', 1, 'MarkerSize', 8);
plot(tvec*1e9, results_mex.wy*1e6, 'c.-', ...
    'LineWidth', 1, 'MarkerSize', 8);
plot(tvec*1e9, results_class.wy*1e6, 'm.--', ...
    'LineWidth', 1, 'MarkerSize', 8);
xlabel('Time (ns)'); ylabel('w_{1/e^2} (\mum)');
title('Beam width vs time');
legend('MEX w_x', 'Class w_x', 'MEX w_y', 'Class w_y', 'Location', 'best');
grid on;

% --- 5. Radius of curvature vs time ---
subplot(2, 3, 5);
plot(tvec*1e9, results_mex.Rx*1e3, 'b.-', ...
    'LineWidth', 1, 'MarkerSize', 8); hold on;
plot(tvec*1e9, results_class.Rx*1e3, 'r.--', ...
    'LineWidth', 1, 'MarkerSize', 8);
plot(tvec*1e9, results_mex.Ry*1e3, 'c.-', ...
    'LineWidth', 1, 'MarkerSize', 8);
plot(tvec*1e9, results_class.Ry*1e3, 'm.--', ...
    'LineWidth', 1, 'MarkerSize', 8);
xlabel('Time (ns)'); ylabel('R (mm)');
title('Radius of curvature vs time');
legend('MEX R_x', 'Class R_x', 'MEX R_y', 'Class R_y', 'Location', 'best');
grid on;

% --- 6. Power vs time ---
% Power = integral of instantaneous intensity over the transverse plane.
P_t = 0.5 * epsilon_0 * c * squeeze(sum(sum(abs(E_pulsed).^2, 1), 2)) * dx * dy;

subplot(2, 3, 6);
plot(tvec*1e9, P_t*1e-3, 'b.-', 'LineWidth', 1.5, 'MarkerSize', 8);
xlabel('Time (ns)'); ylabel('Power (kW)');
title(sprintf('Power vs time\nPeak = %.2f kW', max(P_t)*1e-3));
grid on;

sgtitle(sprintf('Pulsed Gaussian beam:  FWHM_{beam}=%.3f mm,  FWHM_{pulse}=%.1f ns,  %.0f mJ', ...
    FWHM*1e3, FWHM_t*1e9, pulse_energy*1e3), 'FontWeight', 'bold');
