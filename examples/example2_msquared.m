%% Example2: M^2 measurement of a Gaussian beam using msquared_mex and Msquared
%
% Same beam construction as example_msquared.m, but computes the beam
% parameters with both the MEX function and the native MATLAB class and
% prints the two sets of numbers side-by-side for easy comparison.
%
% Grid:    3 mm x 3 mm, 128 x 128 points
% Beam:    FWHM diameter = 0.750 mm,  R = 300 mm
% Wavelength: 1064 nm (Nd:YAG)

%% Physical parameters (all lengths in mm)
epsilon_0 = 8.85e-12;           % vacuum permittivity
c       = 3e8;                  % speed of light
Lx      = 3.0e-3;               % grid side length [m]
Nx      = 128;                  % points per side (total N x N)
dx      = Lx / (Nx-1);          % x step [mm]
Ly      = Lx;
Ny      = Nx;
dy      = dx;

FWHM    = 0.750e-3;             % beam full-width at half-max diameter [mm]
R       = 300e-3;               % radius of curvature [mm]
lambda  = 1064e-9;              % wavelength [m]  (1064 nm)

k0      = 2 * pi / lambda;      % wavenumber [rad/mm]

fprintf('Grid:      %d x %d points over %.2f x %.2f mm\n', Nx, Nx, Lx, Lx);
fprintf('dx:        %.3f um\n', dx*1e6);
fprintf('FWHM (dia) %.3f mm\n', FWHM*1e3);
fprintf('roc:       %.3f mm\n', R*1e3);
fprintf('lambda:    %.0f nm\n', lambda*1e9);

%% Build grid (MATLAB meshgrid convention: X varies along columns, Y varies along rows)
%  msquared_mex expects: xgrid varies in the 1st dim (rows), ygrid varies in
%  the 2nd dim (columns).  We can pass 1-D vectors and let the MEX build the
%  meshgrid internally with the correct orientation, or generate the 2
%  meshgrid arrays. We'll do the second:
xvec    = linspace(-Lx/2, Lx/2, Nx);  % vector of x positions, centered about 0
yvec    = xvec; %
[ymat, xmat] = meshgrid(yvec, xvec);

%% Construct Gaussian field:  E(x,y) = exp(-r^2/w^2) * exp(-i*k*r^2/(2R))
%
%  w : 1/e^2 intensity radius (field amplitude falls as exp(-r^2/w^2))
%      FWHM = w * sqrt(2*ln2)  -->  w = FWHM / sqrt(2*ln2)

w = FWHM / sqrt(2 * log(2));            % 1/e^2 radius [mm]
rSq = xmat.^2 + ymat.^2;                % squared radial distance [mm^2]

% Amplitude factor: Gaussian envelope exp(-r^2/w^2)
amplitude = exp(-rSq / w^2);

% Phase factor: spherical wavefront exp(-i*k*r^2/(2R))
phase = exp(-1i * k0 * rSq / (2 * R));

% Complex field
E = amplitude .* phase;

fprintf('\nBeam 1/e^2 radius w = %.4f mm\n', w*1e3);

%% Normalise field power to 1 W (optional, but good practice)
totalPower = 0.5*epsilon_0*c*trapz(trapz(abs(E).^2)) * dx*dy;
E = E / sqrt(totalPower);   % now sum(|E|^2)*dx^2 = 1 W
fprintf('Total power (after normalisation) = %.2f W\n', 0.5*epsilon_0*c*trapz(trapz(abs(E).^2)) * dx*dy);

%% Call msquared_mex
%  Input: field (Nx x Ny complex), xgrid (Nx x 1), ygrid (1 x Ny), lambda
%  Output: struct with M^2, beam widths, RoC, etc.

fprintf('\nRunning msquared_mex ...\n');
tic;
results_mex = beam.msquared_mex(E, xmat, ymat, lambda);
t_mex = toc;

%% Call Msquared
fprintf('Running Msquared.calculate() ...\n');
m2 = beam.Msquared();
m2.field      = E;
m2.xgrid      = xmat;
m2.ygrid      = ymat;
m2.wavelength = lambda;
tic;
results_class = m2.calculate();
t_class = toc;

%% Display results side-by-side

fprintf('\n========== Beam analysis results (MEX vs Msquared) ==========\n');
fprintf('Quantity          MEX                MATLAB class       Units\n');
fprintf('%s\n', repmat('-', 1, 70));
fprintf('M^2_x:            %.4f             %.4f\n',    results_mex.M2_x,    results_class.M2_x);
fprintf('M^2_y:            %.4f             %.4f\n',    results_mex.M2_y,    results_class.M2_y);
fprintf('wx (1/e^2):       %.4f mm          %.4f mm\n', results_mex.wx*1e3,  results_class.wx*1e3);
fprintf('wy (1/e^2):       %.4f mm          %.4f mm\n', results_mex.wy*1e3,  results_class.wy*1e3);
fprintf('wx0 (waist):      %.4f mm          %.4f mm\n', results_mex.wx0*1e3, results_class.wx0*1e3);
fprintf('wy0 (waist):      %.4f mm          %.4f mm\n', results_mex.wy0*1e3, results_class.wy0*1e3);
fprintf('Rx:               %.1f mm          %.1f mm\n', results_mex.Rx*1e3,  results_class.Rx*1e3);
fprintf('Ry:               %.1f mm          %.1f mm\n', results_mex.Ry*1e3,  results_class.Ry*1e3);
fprintf('z0x (Rayleigh):   %.2f mm          %.2f mm\n', results_mex.z0x*1e3, results_class.z0x*1e3);
fprintf('z0y (Rayleigh):   %.2f mm          %.2f mm\n', results_mex.z0y*1e3, results_class.z0y*1e3);
fprintf('xBar centroid:    %.6f mm        %.6f mm\n', results_mex.xBar*1e3, results_class.xBar*1e3);
fprintf('yBar centroid:    %.6f mm        %.6f mm\n', results_mex.yBar*1e3, results_class.yBar*1e3);
fprintf('kxBar centroid:   %.6f rad/mm     %.6f rad/mm\n', results_mex.kxBar*1e-3, results_class.kxBar*1e-3);
fprintf('kyBar centroid:   %.6f rad/mm     %.6f rad/mm\n', results_mex.kyBar*1e-3, results_class.kyBar*1e-3);
fprintf('wkx (angular):    %.6f rad/mm     %.6f rad/mm\n', results_mex.wkx*1e-3, results_class.wkx*1e-3);
fprintf('wky (angular):    %.6f rad/mm     %.6f rad/mm\n', results_mex.wky*1e-3, results_class.wky*1e-3);
fprintf('===================================================================\n');
fprintf('MEX time:   %.4f s\n', t_mex);
fprintf('Class time: %.4f s\n', t_class);

%% Visualisation

fh = figure('Position', [100 100 1200 500]);

% --- Intensity profile ---
ax1 = axes(fh);
ax2 = axes(fh);
ax3 = axes(fh);
subplot(1, 3, 1, ax1);
[~,a] = contourf(xmat*1e6, ymat*1e6, 0.5*epsilon_0*c*abs(E).^2);
axis image; axis xy;
xlabel('x (\mum)'); ylabel('y (\mum)');
title('Irradiance [W/m^2]');
colorbar('Location', 'eastoutside');
colormap(gca, 'hot');

% --- Phase profile ---
subplot(1, 3, 2, ax2);
[~,a2] = contourf(xmat*1e6, ymat*1e6, angle(E));
axis image; axis xy;
xlabel('x (\mum)'); ylabel('y (\mum)');
title('Phase (rad)');
colorbar('Location', 'eastoutside');
colormap(gca, 'jet');
clim([-pi, pi]);

% --- Cross-section with Gaussian fit ---
subplot(1, 3, 3, ax3);
mid = Nx/2 + 1;
r_mm = xvec;
I_cut = 0.5*epsilon_0*c*abs(E(:, mid)).^2;          % vertical cut through centre
I_max = max(I_cut);
I_fwhm_half = I_max / 2;

% Find FWHM from the cut
above_half = I_cut >= I_fwhm_half;
% Edge positions
idx_left  = find(above_half, 1, 'first');
idx_right = find(above_half, 1, 'last');
FWHM_measured = (xvec(idx_right) - xvec(idx_left));

plot(r_mm*1e3, I_cut, 'b-', 'LineWidth', 1.5); hold on;
yline(I_fwhm_half, 'r--', 'LineWidth', 1);
xline([xvec(idx_left), xvec(idx_right)]*1e3, 'g--', 'LineWidth', 1);
xlabel('Position (\mum)'); ylabel('Intensity');
title(sprintf('Cross-section\nFWHM = %.3f mm (requested %.3f mm)', ...
    FWHM_measured, FWHM));
legend('Intensity', 'Half-max', 'FWHM edges', 'Location', 'north');
grid on;

%% Check: for an ideal Gaussian, M^2 = 1 and measured R should match

fprintf('\nIdeal Gaussian check:\n');
fprintf('  M^2 should be ~1:  M2_x = %.4f / %.4f,  M2_y = %.4f / %.4f\n', ...
    results_mex.M2_x, results_class.M2_x, results_mex.M2_y, results_class.M2_y);
fprintf('  Measured R should be ~%.1f mm:  Rx = %.1f / %.1f mm,  Ry = %.1f / %.1f mm\n', ...
    R*1e3, results_mex.Rx*1e3, results_class.Rx*1e3, results_mex.Ry*1e3, results_class.Ry*1e3);
fprintf('  Measured w should be ~%.4f mm:  wx = %.4f / %.4f mm,  wy = %.4f / %.4f mm\n', ...
    w*1e3, results_mex.wx*1e3, results_class.wx*1e3, results_mex.wy*1e3, results_class.wy*1e3);
