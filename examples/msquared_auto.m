function results = msquared_auto(field, xgrid, ygrid, wavelength, varargin)
%MSQUARED_AUTO Auto-detect compiled MEX or fall back to MATLAB class.
%   RESULTS = MSQUARED_AUTO(FIELD, XGRID, YGRID, WAVELENGTH)
%   RESULTS = MSQUARED_AUTO(..., 'pulse_flat')
%   RESULTS = MSQUARED_AUTO(..., 'second_moment')
%   RESULTS = MSQUARED_AUTO(..., 'pulse')
%
%   If beam.msquared_mex is compiled and on the path, calls it directly.
%   Otherwise, falls back to beam.Msquared() for full compatibility.

    % Use which() instead of exist(...,'file') because MATLAB's exist
    % returns 0 for namespaced MEX files (e.g. +beam/msquared_mex.mexw64).
    if ~isempty(which('beam.msquared_mex'))
        results = beam.msquared_mex(field, xgrid, ygrid, wavelength, varargin{:});
        return;
    end

    m2 = beam.Msquared();
    m2.field = field;
    m2.xgrid = xgrid;
    m2.ygrid = ygrid;
    m2.wavelength = wavelength;

    if nargin > 4 && ischar(varargin{1})
        mode = varargin{1};
        switch mode
            case {'pulse_flat', 'calculate_pulse_flat'}
                results = m2.calculate_pulse_flat();
            case 'second_moment'
                results = m2.calculate_second_moment();
            case 'pulse'
                results = m2.calculate();   % no-op alias
            otherwise
                error('Unknown mode: %s', mode);
        end
    else
        results = m2.calculate();
    end
end
