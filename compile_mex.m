%% compile_mex.m
%
% Compile +beam/private/msquared_mex.cpp with
% platform-specific optimization flags.  Run this from the repo root.
%
% Example:
%   >> compile_mex
%
% Platform auto-detection:
%   - Windows (MSVC or MinGW): compiles directly in MATLAB with
%     optimized C++ flags and links MATLAB's FFTW library.
%   - Linux / WSL: compiles directly in MATLAB with g++,
%     system FFTW3, OpenMP, and native AVX2.
%
% The existing binaries are backed up before overwriting.

function compile_mex(profile)
    if nargin < 1
        profile = false;
    end

    %% macOS: build via system mex (Apple Clang, MATLAB's bundled libomp
    %  + libmwfftw3.dylib).
    %
    %  Mirrors the Linux/WSL block: uses shell-mode mex so CXXOPTIMFLAGS
    %  is passed through reliably.  The configured C++ compiler is whatever
    %  `mex -setup C++` selected (Apple Clang by default on R2023b+).
    %
    %  OpenMP notes:
    %   * Apple Clang 16+ removed the `-fopenmp` driver flag.  We pass the
    %     LLVM frontend flag directly via `-Xclang -fopenmp`.
    %   * MATLAB R2023b+ ships its own libomp.dylib and omp.h; we link
    %     those (no Homebrew / Xcode CLT libomp dependency).  If a
    %     particular MATLAB install lacks the bundled libomp, we fall back
    %     to a no-OpenMP build with a warning.
    %
    %  FFTW: linked directly from MATLAB's runtime dylib; -Wl,-rpath
    %  embeds the search path so the .mexmaca64 / .mexmaci64 loads without
    %  DYLD_LIBRARY_PATH gymnastics at run time.
    %
    %  Detects both Apple Silicon (maca64) and Intel (maci64) MATLAB.
    if ismac
        repoDir = fileparts(mfilename('fullpath'));
        outputDir = fullfile(repoDir, '+beam');
        if ~exist(outputDir, 'dir')
            mkdir(outputDir);
        end

        arch = computer('arch');   % 'maca64' on Apple Silicon, 'maci64' on Intel

        % --- Locate MATLAB's bundled libomp (optional) and libmwfftw3 (required) ---
        ompInclude = fullfile(matlabroot, 'toolbox', 'eml', ...
                              'externalDependency', 'omp', arch, 'include');
        ompLib     = fullfile(matlabroot, 'bin', arch, 'libomp.dylib');
        fftwLib    = fullfile(matlabroot, 'bin', arch, 'libmwfftw3.dylib');

        if ~exist(fftwLib, 'file')
            error('MATLAB FFTW dylib not found: %s', fftwLib);
        end

        hasOmp = exist(ompLib, 'file') == 2 && exist(ompInclude, 'dir') == 7;
        if ~hasOmp
            warning('compile_mex:noOMP', ...
                'MATLAB libomp not found at %s; building without OpenMP.', ompLib);
        end

        % Both libomp (when present) and libmwfftw3 live under
        % <matlabroot>/bin/<arch>, so a single rpath covers both.
        rpathArg = sprintf('-Wl,-rpath,%s', fileparts(fftwLib));

        % Clang on macOS accepts the same -ffast-math family as GCC;
        % -march=native resolves to the host armv8.x variant (or x86_64 on
        % Intel Macs).  -funroll-loops is a no-op for Clang (default) but
        % kept for symmetry with the Linux build.
        ccFlags = '-O3 -march=native -ffast-math -funroll-loops -DNDEBUG';

        if hasOmp
            ccFlags = [ccFlags, ' -Xclang -fopenmp'];
        end

        if profile
            ccFlags = [ccFlags, ' -DMEX_PROFILE'];
            fprintf('Profiling enabled.\n');
        end

        % On macOS, a user-set LDFLAGS env var REPLACES the mexopts' default
        % LDFLAGS (unlike Linux, where it merges).  The mexopts' default
        % supplies -bundle (critical: makes ld produce a Mach-O bundle
        % rather than an executable expecting _main) and -stdlib=libc++
        % (the source is C++), so we re-add them explicitly here.
        if hasOmp
            ompIncludeArg = ['-I' ompInclude];
            ompLinkArg    = ['-bundle -stdlib=libc++ -L' fileparts(ompLib) ' -lomp'];
            cxxFlagsEnv   = ['CXXFLAGS="$CXXFLAGS ' ompIncludeArg '" '];
            ldFlagsEnv    = ['LDFLAGS="$LDFLAGS ' ompLinkArg ' ' rpathArg '" '];
        else
            cxxFlagsEnv   = '';
            ldFlagsEnv    = ['LDFLAGS="$LDFLAGS -bundle -stdlib=libc++ ' rpathArg '" '];
        end

        % Use the absolute path to mex — on macOS MATLAB is normally launched
        % from the GUI and isn't on the user's shell PATH, so a bare `mex`
        % here would fail with "command not found" under zsh/bash.
        mexBin = fullfile(matlabroot, 'bin', arch, 'mex');
        cmd = ['cd "' repoDir '" && "' mexBin '" ' ...
            cxxFlagsEnv ...
            'CXXOPTIMFLAGS="' ccFlags '" ' ...
            'COPTIMFLAGS="'   ccFlags '" ' ...
            ldFlagsEnv ...
            '-outdir "+beam" ' ...
            '-output msquared_mex ' ...
            '+beam/private/msquared_mex.cpp ' ...
            fftwLib];

        if hasOmp
            fprintf('macOS detected — building with mex (%s, libmwfftw3 + libomp)...\n', arch);
        else
            fprintf('macOS detected — building with mex (%s, libmwfftw3, no OpenMP)...\n', arch);
        end
        fprintf('C/C++ optimization flags: %s\n', ccFlags);
        fprintf('FFTW library: %s\n', fftwLib);
        [status, result] = system(cmd);
        fprintf('%s', result);
        if status ~= 0
            error('mex exited with code %d', status);
        end
        fprintf('Build complete.\n');
        return;
    end

    %% Linux / WSL: build via system mex (shell-mode respects CXXOPTIMFLAGS)
    %
    % Default: no direct FFTW; uses mexCallMATLAB("fft2") which gets MATLAB's
    % MKL.  On WSL2 this is ~1.5x, while direct FFTW is also ~1.5x but adds a
    % library dependency, so MKL is the simpler default.
    %
    % To benchmark direct FFTW, uncomment the USE_FFTW block below.
    if isunix && ~ismac
        repoDir = fileparts(mfilename('fullpath'));
        outputDir = fullfile(repoDir, '+beam');
        if ~exist(outputDir, 'dir')
            mkdir(outputDir);
        end

        % --- Default MKL build ---
        ccFlags = '-O3 -march=native -ffast-math -funroll-loops -fopenmp -DNDEBUG';

        % --- Optional FFTW benchmark build ---
        % ccFlags = '-O3 -march=native -ffast-math -funroll-loops -fopenmp -DNDEBUG -DUSE_FFTW';
        % fftwLib = fullfile(matlabroot, 'bin', 'glnxa64', 'libmwfftw3.so');

        if profile
            ccFlags = [ccFlags, ' -DMEX_PROFILE'];
            fprintf('Profiling enabled.\n');
        end

        % Shell-mode mex is the only reliable way to override CXXOPTIMFLAGS on Linux
        cmd = ['cd "' repoDir '" && mex ' ...
            'CXXOPTIMFLAGS="' ccFlags '" ' ...
            'COPTIMFLAGS="' ccFlags '" ' ...
            'LDFLAGS="$LDFLAGS -fopenmp" ' ...
            '-outdir "+beam" ' ...
            '-output msquared_mex ' ...
            '+beam/private/msquared_mex.cpp'];

        % If FFTW is uncommented, also append LINKLIBS to cmd above:
        % cmd = [cmd, ' LINKLIBS="$LINKLIBS -l:' fftwLib '" '];

        fprintf('Linux/WSL detected — building with mex (g++, MKL fft2 + OpenMP)...\n');
        fprintf('C/C++ optimization flags: %s\n', ccFlags);
        [status, result] = system(cmd);
        fprintf('%s', result);
        if status ~= 0
            error('mex exited with code %d', status);
        end
        fprintf('Build complete.\n');
        return;
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
    %  +beam/private/msquared_mex.cpp uses the separate real/imaginary API (R2017b).
    sources = {'+beam/private/msquared_mex.cpp'};
    extraFlags = {{}};   % per-source additional mex arguments

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

    %% MATLAB FFTW library (for direct FFT in +beam/private/msquared_mex.cpp)
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
