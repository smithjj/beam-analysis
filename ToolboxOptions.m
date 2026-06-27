function opts = ToolboxOptions()
%TOOLBOXOPTIONS Packaging options for the Msquared Toolbox.
%   OPTS = ToolboxOptions() returns a matlab.addons.toolbox.ToolboxOptions
%   object configured to package the M-squared beam-quality toolbox.
%
%   Use:
%       opts = ToolboxOptions();
%       matlab.addons.toolbox.packageToolbox(opts);
%
%   See also matlab.addons.toolbox.packageToolbox.

    root = fileparts(mfilename('fullpath'));
    opts = matlab.addons.toolbox.ToolboxOptions(root, ...
        "Identifier", "4fc9c860-046b-492a-868e-27a83a9f6047");

    opts.ToolboxName        = "Msquared Toolbox";
    opts.ToolboxVersion     = "1.0.0";
    opts.AuthorName         = "Jesse Smith";
    opts.AuthorEmail        = "jesse.smith@as-photonics.com";
    opts.AuthorCompany      = "AS-Photonics";
    opts.Summary            = "M-squared beam quality calculator with MEX accelerators";
    opts.Description        = ["Calculate the M-squared beam quality factor and related " ...
                               "beam quantities for 2-D and 3-D electric fields. " ...
                               "Includes a reference MATLAB class (beam.Msquared) and " ...
                               "compiled MEX accelerators (beam.msquared_mex, " ...
                               "beam.msquared_calculate_mex)."];
    opts.MinimumMatlabRelease = "R2020b";
    opts.Platforms          = [];
    opts.SupportedOperations  = ...
        struct(Install = true, Download = true, Update = true, Query = true);

    % Include the top-level README as toolbox documentation.
    % MATLAB will pick up Contents.m and GettingStarted.m automatically.
    opts.ToolboxGettingStartedFile = fullfile(root, "GettingStarted.m");

    % Package output path
    opts.OutputFile = fullfile(root, "MsquaredToolbox.mltbx");
end
