%% Example: Pulsed M^2 measurement of a Gaussian beam using msquared_mex
%
% Constructs a Gaussian transverse field and multiplies by a Gaussian
% temporal pulse envelope, then calls msquared_mex in 'pulse' mode to
% compute per-time-slice beam propagation parameters.
%
% Grid:       3 mm x 3 mm, 128 x 128 points
% Beam:       FWHM diameter = 0.750 mm,  R = 300 mm
% Wavelength: 1064 nm
% Pulse:      FWHM = 3 ns,  t = -10 .. +10 ns, 32 pts, 1 mJ energy

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

%% Call msquared_mex in pulse mode
%
%  Input:  field (Nx x Ny x Nt), xgrid, ygrid, wavelength, 'pulse'
%  Output: struct with standard field names (M2_x, wx, Rx, ...) — each is an Nt x 1 array
%          ('pulse' is accepted as a no-op alias for the default per-slice mode.)
%
fprintf('\nRunning msquared_mex (instantaneous mode) ...\n');
tic;
results = msquared_mex(E_pulsed, xv, yv, lambda, 'pulse');
toc;

%% --- Display results ---

% Centre time slice (closest to t = 0)
[~, t0_idx] = min(abs(tvec));

fprintf('\n========== Beam analysis @ t = %.2f ns ==========\n', tvec(t0_idx)*1e9);
fprintf('M^2_x:          %.4f\n',    results.M2_x(t0_idx)   );
fprintf('M^2_y:          %.4f\n',    results.M2_y(t0_idx)   );
fprintf('wx (1/e^2):     %.4f mm\n', results.wx(t0_idx)*1e3 );
fprintf('wy (1/e^2):     %.4f mm\n', results.wy(t0_idx)*1e3 );
fprintf('wx0 (waist):    %.4f mm\n', results.wx0(t0_idx)*1e3);
fprintf('wy0 (waist):    %.4f mm\n', results.wy0(t0_idx)*1e3);
fprintf('Rx:             %+.1f mm\n',results.Rx(t0_idx)*1e3 );
fprintf('Ry:             %+.1f mm\n',results.Ry(t0_idx)*1e3 );
fprintf('z0x (Rayleigh): %.2f mm\n', results.z0x(t0_idx)*1e3);
fprintf('z0y (Rayleigh): %.2f mm\n', results.z0y(t0_idx)*1e3);
fprintf('=================================================\n');

%% Fluence-based pulse analysis: calculate_pulse_flat equivalent
%  This integrates over time first and returns one scalar result for the
%  whole pulse, rather than one result per time slice.

fprintf('\nRunning msquared_mex (fluence-based pulse_flat mode) ...\n');
tic;
results_flat = msquared_mex(E_pulsed, X2d, Y2d, lambda, 'pulse_flat');
toc;

fprintf('\n Fluence-based pulse analysis %s', newline);
fprintf('pulse_M2_x:      %.4f\n',      results_flat.pulse_M2_x);
fprintf('pulse_M2_y:      %.4f\n',      results_flat.pulse_M2_y);
fprintf('pulse_wx:        %.4f mm\n',   results_flat.pulse_wx*1e3);
fprintf('pulse_wy:        %.4f mm\n',   results_flat.pulse_wy*1e3);
fprintf('pulse_Rx:        %+.1f mm\n',  results_flat.pulse_Rx*1e3);
fprintf('pulse_Ry:        %+.1f mm\n',  results_flat.pulse_Ry*1e3);
fprintf('==================================================%s', newline);

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
plot(tvec*1e9, results.M2_x, 'b.-', ...
    'LineWidth', 1, 'MarkerSize', 8); hold on;
plot(tvec*1e9, results.M2_y, 'r.-', ...
    'LineWidth', 1, 'MarkerSize', 8);
xlabel('Time (ns)'); ylabel('M^2');
title('M^2 vs time');
legend('M^2_x', 'M^2_y', 'Location', 'best');
ylim([min(results.M2_x)*0.999, max(results.M2_x)*1.001]);
grid on;

% --- 4. Beam width vs time ---
subplot(2, 3, 4);
plot(tvec*1e9, results.wx*1e6, 'b.-', ...
    'LineWidth', 1, 'MarkerSize', 8); hold on;
plot(tvec*1e9, results.wy*1e6, 'r.-', ...
    'LineWidth', 1, 'MarkerSize', 8);
xlabel('Time (ns)'); ylabel('w_{1/e^2} (\mum)');
title('Beam width vs time');
legend('w_x', 'w_y', 'Location', 'best');
grid on;

% --- 5. Radius of curvature vs time ---
subplot(2, 3, 5);
plot(tvec*1e9, results.Rx*1e3, 'b.-', ...
    'LineWidth', 1, 'MarkerSize', 8); hold on;
plot(tvec*1e9, results.Ry*1e3, 'r.-', ...
    'LineWidth', 1, 'MarkerSize', 8);
xlabel('Time (ns)'); ylabel('R (mm)');
title('Radius of curvature vs time');
legend('R_x', 'R_y', 'Location', 'best');
grid on;

% --- 6. Check energy in flattened (curvature-removed) field ---
E_flat = results.flattened_Exyz;       % Nx*Ny x Nt, column-major
E_flat = reshape(E_flat, Nx, Nx, Nt);        % Nx x Ny x Nt
energy_flat = 0.5 * epsilon_0 * c * sum(abs(E_flat(:)).^2) * dx * dy * dt;

I_slice = squeeze(sum(sum(abs(E_flat).^2, 1), 2));
I_slice = I_slice / max(I_slice);
subplot(2, 3, 6);
plot(tvec*1e9, I_slice, 'k.-', 'LineWidth', 1.5, 'MarkerSize', 8);
xlabel('Time (ns)'); ylabel('Norm. energy per slice');
title(sprintf('Flat-field envelope\nTotal = %.3f mJ', energy_flat*1e3));
grid on;

sgtitle(sprintf('Pulsed Gaussian beam:  FWHM_{beam}=%.3f mm,  FWHM_{pulse}=%.1f ns,  %.0f mJ', ...
    FWHM*1e3, FWHM_t*1e9, pulse_energy*1e3), 'FontWeight', 'bold');
