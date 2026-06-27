%% compile_mex.m
%
% Compile msquared_mex.cpp and msquared_calculate_mex.cpp with
% platform-specific optimization flags.  Run this from the repo root.
%
% Example:
%   >> compile_mex
%
% The script detects the active C++ compiler and selects appropriate flags:
%   - Microsoft Visual C++ 2022: /O2 /openmp /arch:AVX2 /fp:fast /DNDEBUG
%   - MinGW-w64:                 -O3 -march=native -fopenmp -ffast-math -DNDEBUG
%
% The existing .mexw64 files are backed up before overwriting.

function compile_mex(profile)
    if nargin < 1
        profile = false;
    end

    outputDir = '+beam';
    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end

    %% Detect active compiler
    cfg = mex.getCompilerConfigurations('C++', 'Selected');
    if isempty(cfg)
        error('No C++ MEX compiler is selected. Run ''mex -setup C++'' first.');
    end
    compilerName = cfg.Name;
    fprintf('Selected compiler: %s\n', compilerName);

    %% Choose optimization flags
    if contains(compilerName, 'Microsoft Visual C++', 'IgnoreCase', true)
        % MSVC flags
        ccFlags  = '/O2 /openmp /arch:AVX2 /fp:fast /DNDEBUG';
        ldFlags  = '/O2 /DNDEBUG';
        isMSVC   = true;
    elseif contains(compilerName, 'MinGW', 'IgnoreCase', true)
        % MinGW-w64 flags
        ccFlags  = '-O3 -march=native -fopenmp -ffast-math -DNDEBUG';
        ldFlags  = '-O3 -fopenmp -DNDEBUG';
        isMSVC   = false;
    else
        warning('Unknown compiler ''%s''; using generic -O2 flags.', compilerName);
        ccFlags  = '-O2 -DNDEBUG';
        ldFlags  = '-O2 -DNDEBUG';
        isMSVC   = false;
    end

    if profile
        ccFlags = [ccFlags, ' -DMEX_PROFILE'];
        fprintf('Profiling enabled.\n');
    end

    fprintf('C/C++ optimization flags: %s\n', ccFlags);
    fprintf('Linker optimization flags: %s\n', ldFlags);

    %% Source files to compile
    %  msquared_calculate_mex.cpp uses mxComplexDouble/mxGetComplexDoubles,
    %  which require the R2018a interleaved-complex API.  msquared_mex.cpp
    %  still uses the separate real/imaginary API (R2017b).
    sources = {'msquared_mex.cpp', 'msquared_calculate_mex.cpp'};
    extraFlags = {{}, {'-R2018a'}};   % per-source additional mex arguments

    %% Backup (and temporarily remove) existing .mexw64 files
    %  On Windows a loaded/locked MEX cannot be overwritten.  We move it
    %  aside so the linker can create a fresh file; if compilation fails
    %  we restore the moved file.
    backupSuffix = datestr(now, 'yyyymmdd_HHMMSS');
    origNames  = {};
    movedNames = {};
    for i = 1:numel(sources)
        [~, base, ~] = fileparts(sources{i});
        mexFile = fullfile(outputDir, [base, '.', mexext]);
        if exist(mexFile, 'file')
            backupName = sprintf('%s_backup_%s.%s', base, backupSuffix, mexext);
            copyfile(mexFile, backupName);
            fprintf('Backed up %s -> %s\n', mexFile, backupName);
            % Move original out of the way so linker can write a new file
            tempMoved = sprintf('%s_backup_moved_%s.%s', base, backupSuffix, mexext);
            movefile(mexFile, tempMoved);
            origNames{end+1}  = mexFile;  %#ok<AGROW>
            movedNames{end+1} = tempMoved; %#ok<AGROW>
            fprintf('Moved %s -> %s (will restore on failure)\n', mexFile, tempMoved);
        end
    end

    %% MATLAB FFTW library (for direct FFT in msquared_mex.cpp)
    fftwLib = fullfile(matlabroot, 'lib', 'win64', 'libmwfftw3.lib');
    if ~exist(fftwLib, 'file')
        error('MATLAB FFTW import library not found: %s', fftwLib);
    end

    %% Compile each source file
    for i = 1:numel(sources)
        src = sources{i};
        flags = extraFlags{i};
        [~, base, ~] = fileparts(src);
        useFFTW = strcmp(base, 'msquared_mex');
        thisCcFlags = ccFlags;
        if useFFTW
            thisCcFlags = [thisCcFlags, ' -DUSE_FFTW'];
        end
        fprintf('\nCompiling %s ...\n', src);
        try
            if isMSVC
                if useFFTW
                    mex('-v', '-outdir', outputDir, src, flags{:}, fftwLib, ...
                        ['CXXOPTIMFLAGS=', thisCcFlags], ...
                        ['COPTIMFLAGS=',   thisCcFlags], ...
                        ['LDOPTIMFLAGS=',  ldFlags]);
                else
                    mex('-v', '-outdir', outputDir, src, flags{:}, ...
                        ['CXXOPTIMFLAGS=', thisCcFlags], ...
                        ['COPTIMFLAGS=',   thisCcFlags], ...
                        ['LDOPTIMFLAGS=',  ldFlags]);
                end
            else
                if useFFTW
                    mex('-v', '-outdir', outputDir, src, flags{:}, fftwLib, ...
                        ['CXXOPTIMFLAGS=', thisCcFlags], ...
                        ['COPTIMFLAGS=',   thisCcFlags], ...
                        ['LDFLAGS=',       ldFlags]);
                else
                    mex('-v', '-outdir', outputDir, src, flags{:}, ...
                        ['CXXOPTIMFLAGS=', thisCcFlags], ...
                        ['COPTIMFLAGS=',   thisCcFlags], ...
                        ['LDFLAGS=',       ldFlags]);
                end
            end
            fprintf('  -> %s built successfully.\n', mexFile);
        catch ME
            fprintf('  -> ERROR building %s:\n%s\n', src, ME.message);
            % Restore any moved originals so the repo is left in a usable state
            for j = 1:numel(origNames)
                if exist(movedNames{j}, 'file') && ~exist(origNames{j}, 'file')
                    movefile(movedNames{j}, origNames{j});
                    fprintf('Restored %s -> %s\n', movedNames{j}, origNames{j});
                end
            end
            rethrow(ME);
        end
    end

    % Clean up temporary moved files (keep timestamped backups)
    for j = 1:numel(movedNames)
        if exist(movedNames{j}, 'file')
            delete(movedNames{j});
        end
    end

    fprintf('\nAll MEX files compiled.\n');
end
