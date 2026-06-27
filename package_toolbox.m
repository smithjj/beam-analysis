% package_toolbox.m
% Package the Msquared Toolbox as an installable .mltbx file.
%
% Run from the repository root in MATLAB:
%   >> package_toolbox
%
% See also ToolboxOptions, matlab.addons.toolbox.packageToolbox.

opts = ToolboxOptions();
matlab.addons.toolbox.packageToolbox(opts);
