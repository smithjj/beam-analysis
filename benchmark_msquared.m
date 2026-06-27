%% benchmark_msquared.m
%
% Quick scaling benchmark for msquared_mex vs Msquared_class.
% Builds Gaussian beams of varying grid size and slice count and reports
% timings for the standard per-slice mode and pulse_flat mode.

function benchmark_msquared()
    epsilon_0 = 8.854187817e-12;
    c         = 299792458;

    FWHM     = 0.750e-3;   % m
    R        = 300e-3;     % m
    lambda   = 1064e-9;    % m
    k0       = 2 * pi / lambda;
    w        = FWHM / sqrt(2 * log(2));

    grid_sizes = [64, 128, 256];
    nt_values  = [1, 32, 129];

    fprintf('%-6s %-6s %-14s %-14s %-14s %-14s\n', ...
        'Nx', 'Nt', 'MEX standard', 'Class standard', 'MEX flat', 'Class flat');
    fprintf('%s\n', repmat('-', 1, 80));

    for Nx = grid_sizes
        xv = linspace(-1.5e-3, 1.5e-3, Nx);
        [X, Y] = meshgrid(xv, xv);
        rSq = X.^2 + Y.^2;
        E_cw = exp(-rSq / w^2) .* exp(-1i * k0 * rSq / (2 * R));

        for Nt = nt_values
            %% Build pulsed field
            tvec = linspace(-10e-9, 10e-9, Nt);
            FWHM_t = 3e-9;
            pulse_env = exp(-2 * log(2) * (tvec / FWHM_t).^2);
            E_pulsed = E_cw .* reshape(pulse_env, 1, 1, Nt);

            % Normalize roughly to 1 mJ (not critical for M2)
            dx = xv(2) - xv(1);
            if Nt > 1
                dt = tvec(2) - tvec(1);
            else
                dt = 1;  % arbitrary for single-slice normalization
            end
            energy = 0.5 * epsilon_0 * c * sum(abs(E_pulsed(:)).^2) * dx^2 * dt;
            if energy > 0
                E_pulsed = E_pulsed * sqrt(1e-3 / energy);
            end

            %% Standard per-slice mode
            t_mex = timeit(@() msquared_mex(E_pulsed, xv, xv, lambda, 'pulse'), 1);

            m2 = Msquared_class();
            m2.field = E_pulsed;
            m2.xgrid = xv;
            m2.ygrid = xv;
            m2.wavelength = lambda;
            t_class = timeit(@() m2.calculate(), 1);

            %% pulse_flat mode (only meaningful for Nt > 1)
            if Nt > 1
                t_mex_flat = timeit(@() msquared_mex(E_pulsed, X, Y, lambda, 'pulse_flat'), 1);
                t_class_flat = timeit(@() m2.calculate_pulse_flat(), 1);
            else
                t_mex_flat = NaN;
                t_class_flat = NaN;
            end

            fprintf('%-6d %-6d %-14.4f %-14.4f %-14.4f %-14.4f\n', ...
                Nx, Nt, t_mex, t_class, t_mex_flat, t_class_flat);
        end
        fprintf('\n');
    end
end
