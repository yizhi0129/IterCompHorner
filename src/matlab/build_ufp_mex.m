function build_ufp_mex()
%BUILD_UFP_MEX Compile the binary64 exponent-based ufp implementation.
%   No external libraries are required. Output is placed in src/mex.
    source_dir = fileparts(fileparts(mfilename('fullpath')));
    core_dir = fullfile(source_dir, 'c');
    mex_dir = fullfile(source_dir, 'mex');
    mex('-R2018a', '-O', '-outdir', mex_dir, ['-I', core_dir], ...
        fullfile(mex_dir, 'ufp_mex.c'), fullfile(core_dir, 'ufp.c'));
end
