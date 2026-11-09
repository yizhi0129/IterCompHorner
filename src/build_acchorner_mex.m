function build_acchorner_mex()
%BUILD_ACCHORNER_MEX Compile AccHorner_mex with MPFR/GMP.

    root = fileparts(mfilename('fullpath'));
    source = fullfile(root, 'AccHorner_mex.c');

    mpfrRoot = getenv('MPFR_ROOT');
    if strlength(mpfrRoot) == 0
        mpfrRoot = '/opt/homebrew';
    end

    includeFlag = ['-I', char(fullfile(mpfrRoot, 'include'))];
    libFlag = ['-L', char(fullfile(mpfrRoot, 'lib'))];

    mex('-R2018a', includeFlag, libFlag, '-lmpfr', '-lgmp', source);
end
