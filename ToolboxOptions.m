function opts = ToolboxOptions()
%TOOLBOXOPTIONS Packaging options for the Msquared Toolbox.
%   OPTS = ToolboxOptions() returns a matlab.addons.toolbox.ToolboxOptions
%   object configured to package the M-squared beam-quality toolbox.
%
%   Ship model:
%       The .mltbx ships the MATLAB source (+beam/Msquared.m), the MEX C++
%       source (+beam/private/msquared_mex.cpp), the build script
%       (compile_mex.m), the examples, and the docs.  It does NOT ship
%       prebuilt MEX binaries: those are excluded via opts.ToolboxFiles
%       because they are CPU- and OS-specific.  End users run compile_mex
%       after install to build a binary for their own machine, which lands
%       in +beam/ alongside the shipped source.  This keeps the toolbox
%       itself cross-platform.
%
%   Use:
%       opts = ToolboxOptions();
%       matlab.addons.toolbox.packageToolbox(opts);
%
%   See also matlab.addons.toolbox.packageToolbox, compile_mex, localListToolboxFiles.

    root = fileparts(mfilename('fullpath'));
    opts = matlab.addons.toolbox.ToolboxOptions(root, ...
        "4fc9c860-046b-492a-868e-27a83a9f6047");

    opts.ToolboxName        = "Msquared Toolbox";
    opts.ToolboxVersion     = "1.0.1";
    opts.AuthorName         = "Jesse Smith";
    opts.AuthorEmail        = "jesse.smith@as-photonics.com";
    opts.AuthorCompany      = "AS-Photonics";
    opts.Summary            = 'M-squared beam quality calculator with MEX accelerators';
    opts.Description        = 'Calculate the M-squared beam quality factor and related beam quantities for 2-D and 3-D electric-field arrays. Includes a reference MATLAB class (beam.Msquared) and a compiled MEX accelerator (beam.msquared_mex).';

    % Explicit file list — see localListToolboxFiles below.  Excludes
    % *.mex* (compiled binaries) and the .mltbx artifact itself.
    opts.ToolboxFiles = localListToolboxFiles(root);

    % Package output path
    opts.OutputFile = fullfile(root, "MsquaredToolbox.mltbx");
end

function files = localListToolboxFiles(root)
%LOCALLISTTOOLBOXFILES Enumerate files to ship in the .mltbx.
%   FILES = localListToolboxFiles(ROOT) walks ROOT recursively and returns
%   a cell array of file paths relative to ROOT, excluding:
%     - MEX binaries (*.mexw64, *.mexa64, *.mexmaca64, *.mexmaci64, etc.)
%     - pre-built MEX zip archives (msquared_mex.zip) shipped at the root
%     - the built .mltbx artifact itself
%     - the compile_wsl.sh wrapper (compile_mex.m handles all platforms
%       now, so the shell wrapper is no longer needed in the toolbox)
%     - backup / scratch files left over from compile_mex's move dance
%     - macOS Finder detritus (.DS_Store) and VCS metadata (.git/)

    d = dir(fullfile(root, '**', '*'));
    d = d(~[d.isdir]);

    % Paths relative to the toolbox root.
    absPaths = fullfile({d.folder}, {d.name});
    rootPrefix = [root, filesep];
    n = length(rootPrefix);
    relPaths = cellfun(@(p) p(n+1:end), absPaths, 'UniformOutput', false);

    isMexNonWin = ~cellfun(@isempty, regexp(relPaths, '\.mex(mac|a)\w*$', 'once'));
    isZipMEX    = strcmp(relPaths, 'msquared_mex.zip');
    isMltbx   = strcmp(relPaths, 'MsquaredToolbox.mltbx');
    isShell   = strcmp(relPaths, 'compile_wsl.sh');
    isBackup  = ~cellfun(@isempty, regexp(relPaths, '_backup(_moved)?_', 'once'));
    isDSStore = strcmp(relPaths, '.DS_Store');
    isVCS     = ~cellfun(@isempty, regexp(relPaths, '^\.git(/|$)', 'once'));

    keep = ~(isMexNonWin | isZipMEX | isMltbx | isShell | isBackup | isDSStore | isVCS);
    files = sort(relPaths(keep));
end
