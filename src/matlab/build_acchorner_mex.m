function build_acchorner_mex()
%BUILD_ACCHORNER_MEX Compile AccHorner_mex with MPFR/GMP.

    source_dir = fileparts(fileparts(mfilename('fullpath')));
    core_dir = fullfile(source_dir, 'c');
    mex_dir = fullfile(source_dir, 'mex');

    mpfrRoot = getenv('MPFR_ROOT');
    if strlength(mpfrRoot) == 0
        mpfrRoot = '/opt/homebrew';
    end

    includeFlag = ['-I', char(fullfile(mpfrRoot, 'include'))];
    libFlag = ['-L', char(fullfile(mpfrRoot, 'lib'))];

    mex('-R2018a', '-outdir', mex_dir, ['-I', core_dir], includeFlag, libFlag, ...
        fullfile(mex_dir, 'AccHorner_mex.c'), fullfile(core_dir, 'acchorner.c'), ...
        '-lmpfr', '-lgmp');
end
