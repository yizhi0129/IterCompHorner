function build_ufp_mex()
%BUILD_UFP_MEX Compile the binary64 exponent-based ufp implementation.
%   No external libraries are required. Output is placed beside this file.
    root = fileparts(mfilename('fullpath'));
    mex('-R2018a', '-O', '-outdir', root, fullfile(root, 'ufp_mex.c'));
end
