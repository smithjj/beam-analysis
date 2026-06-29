%% Quick single-run benchmark of inner-OpenMP effect (Nt=1..64)

FWHM    = 0.750e-3;
R       = 300e-3;
lambda  = 1064e-9;
k0      = 2 * pi / lambda;
w       = FWHM / sqrt(2 * log(2));

Nx = 256;
xv = linspace(-1.5e-3, 1.5e-3, Nx);
[X, Y] = meshgrid(xv, xv);
rSq = X.^2 + Y.^2;
E_cw = exp(-rSq / w^2) .* exp(-1i * k0 * rSq / (2 * R));

nt_values = [1, 2, 4, 8, 16, 32, 64];

m2 = beam.Msquared();
m2.xgrid = xv;
m2.ygrid = xv;
m2.wavelength = lambda;

fprintf('%-6s %-14s %-14s %-10s %-10s\n', 'Nt', 'MEX (s)', 'Class (s)', 'Speedup', 'Winner');
fprintf('%s\n', repmat('-', 1, 60));

for Nt = nt_values
    tvec = linspace(-10e-9, 10e-9, Nt);
    FWHM_t = 3e-9;
    pulse_env = exp(-2 * log(2) * (tvec / FWHM_t).^2);
    E_pulsed = E_cw .* reshape(pulse_env, 1, 1, Nt);

    % A few warmup calls
    beam.msquared_mex(E_pulsed, xv, xv, lambda, 'pulse');
    beam.msquared_mex(E_pulsed, xv, xv, lambda, 'pulse');

    nrep = max(1, floor(1 / Nt));  % more reps for small Nt
    tic;
    for k = 1:nrep
        beam.msquared_mex(E_pulsed, xv, xv, lambda, 'pulse');
    end
    t_mex = toc / nrep;

    m2.field = E_pulsed;
    tic;
    for k = 1:nrep
        m2.calculate();
    end
    t_class = toc / nrep;

    speedup = t_class / t_mex;
    winner = t_mex < t_class * 0.95;
    loser  = t_mex > t_class * 1.05;
    if winner
        tag = 'MEX wins';
    elseif loser
        tag = 'Class wins';
    else
        tag = 'Tie';
    end
    fprintf('%-6d %-14.4f %-14.4f %-10.2f %-10s\n', Nt, t_mex, t_class, speedup, tag);
end
