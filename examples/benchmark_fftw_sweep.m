%% benchmark_fftw_sweep
%
% Sweep Nt for fixed Nx=256 and find where MEX-vs-class break-even occurs.
% Helps identify the FFTW/MKL crossover point.

function benchmark_fftw_sweep()
    FWHM    = 0.750e-3;
    R       = 300e-3;
    lambda  = 1064e-9;
    k0      = 2 * pi / lambda;
    w       = FWHM / sqrt(2 * log(2));

    Nx = 64;
    xv = linspace(-1.5e-3, 1.5e-3, Nx);
    [X, Y] = meshgrid(xv, xv);
    rSq = X.^2 + Y.^2;
    E_cw = exp(-rSq / w^2) .* exp(-1i * k0 * rSq / (2 * R));

    % Sweep Nt around the anomaly region
    nt_values = [1, 2, 4, 8, 16, 24, 32, 40, 48, 56, 64, 80, 96, 112, 129, 160, 200];
    fprintf('Nx=Ny=%i\n', Nx);
    fprintf('%-6s %-14s %-14s %-10s %-10s\n', ...
        'Nt', 'MEX (s)', 'Class (s)', 'Speedup', 'Winner');
    fprintf('%s\n', repmat('-', 1, 60));

    for Nt = nt_values
        tvec = linspace(-10e-9, 10e-9, Nt);
        FWHM_t = 3e-9;
        pulse_env = exp(-2 * log(2) * (tvec / FWHM_t).^2);
        E_pulsed = E_cw .* reshape(pulse_env, 1, 1, Nt);

        t_mex = timeit(@() msquared_auto(E_pulsed, xv, xv, lambda, 'pulse'), 1);

        m2 = beam.Msquared();
        m2.field = E_pulsed;
        m2.xgrid = xv;
        m2.ygrid = xv;
        m2.wavelength = lambda;
        t_class = timeit(@() m2.calculate(), 1);

        speedup = t_class / t_mex;
        winner = t_mex < t_class * 0.95;
        loser  = t_mex > t_class * 1.05;
        tag = '';
        if winner
            tag = 'MEX wins';
        elseif loser
            tag = 'Class wins';
        else
            tag = 'Tie';
        end

        fprintf('%-6d %-14.4f %-14.4f %-10.2f %-10s\n', ...
            Nt, t_mex, t_class, speedup, tag);
    end

    fprintf('\n===== pulse_flat (fluence-based) =====\n');
    fprintf('%-6s %-14s %-14s %-10s %-10s\n', ...
        'Nt', 'MEX_flat (s)', 'Class_flat (s)', 'Speedup', 'Winner');
    fprintf('%s\n', repmat('-', 1, 60));

    for Nt = nt_values(2:end)  % skip Nt=1, pulse_flat meaningless
        tvec = linspace(-10e-9, 10e-9, Nt);
        pulse_env = exp(-2 * log(2) * (tvec / FWHM_t).^2);
        E_pulsed = E_cw .* reshape(pulse_env, 1, 1, Nt);

        t_mex = timeit(@() msquared_auto(E_pulsed, X, Y, lambda, 'pulse_flat'), 1);
        m2.field = E_pulsed;
        t_class = timeit(@() m2.calculate_pulse_flat(), 1);

        speedup = t_class / t_mex;
        winner = t_mex < t_class * 0.95;
        loser  = t_mex > t_class * 1.05;
        tag = '';
        if winner
            tag = 'MEX wins';
        elseif loser
            tag = 'Class wins';
        else
            tag = 'Tie';
        end

        fprintf('%-6d %-14.4f %-14.4f %-10.2f %-10s\n', ...
            Nt, t_mex, t_class, speedup, tag);
    end
end
