%% Getting Started with the Msquared Toolbox
% The |Msquared Toolbox| calculates the M-squared beam-quality factor and
% related beam quantities for 2-D (x,y) or 3-D (x,y,t) electric-field
% arrays.  It provides a reference MATLAB class plus compiled MEX
% accelerators.
%
% Add the toolbox root to your MATLAB path, then the |+beam| namespace is
% automatically available:
%
%   addpath('C:\path\to\beam-analysis');
%   m2 = beam.Msquared();          % create calculator object
%   m2.xgrid = xmat;               % set spatial grids
%   m2.ygrid = ymat;
%   m2.field = E;                  % electric field (V/m)
%   results = m2.calculate();      % compute M-squared
%
% For a speed boost, call the MEX replacement directly:
%
%   results = beam.msquared_mex(E, xmat, ymat, lambda);
%
% Open the example scripts in |examples/| for complete use cases, including
% pulsed beams and a performance benchmark.
%
% See also beam.Msquared, beam.msquared_mex, examples.example_msquared.
