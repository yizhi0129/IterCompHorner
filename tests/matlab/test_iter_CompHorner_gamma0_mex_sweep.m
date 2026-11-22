function output_files = test_iter_CompHorner_gamma0_mex_sweep()
%TEST_ITER_COMPHORNER_GAMMA0_MEX_SWEEP Write six relative-error reports.
%   For each degree in [5 10 50 100 300 500], generate one polynomial per
%   target condition 10^(2:35), using rng(exponent) as in the original sweep.
%   Use revised iter_CompHorner_accsum and iter_stop_test_accsum_mex.
%   Report columns and INTLAB reference precision match the original sweep.

    root = fileparts(mfilename('fullpath'));
    if exist('ufp_mex', 'file') ~= 3
        build_ufp_mex();
    end
    saved_rng = rng;
    cleanup_rng = onCleanup(@() rng(saved_rng)); %#ok<NASGU>
    degrees = [5, 10, 50, 100, 300, 500];
    exponents = 2:35;
    output_files = cell(size(degrees));
    for idx = 1:numel(degrees)
        degree = degrees(idx);
        output_files{idx} = fullfile(root, sprintf( ...
            'iter_CompHorner_gamma0_mex_degree%d_sweep_1e%d_1e%d_relerr_output.txt', ...
            degree, exponents(1), exponents(end)));
        run_sweep(degree, exponents, output_files{idx});
    end
end

function run_sweep(degree, exponents, output_file)
    x = 1.03;
    v = 1e-3;
    max_iter = 100;
    fid = fopen(output_file, 'w');
    assert(fid >= 0, 'Cannot open %s.', output_file);
    cleanup_file = onCleanup(@() fclose(fid)); %#ok<NASGU>

    fprintf(fid, '=== iter_CompHorner revised AccSumK-H / MEX-stop condition sweep ===\n\n');
    fprintf(fid, 'implementation = iter_CompHorner_accsum\nstop = iter_stop_test_accsum_mex\n');
    fprintf(fid, 'degree = %d\nx = %.17g\nv = %.17g\nmax_iter = %d\n', ...
        degree, x, v, max_iter);
    fprintf(fid, '%8s  %22s  %9s  %10s  %13s  %14s  %12s\n', ...
        'target', 'actual_condition', 'stopped', 'iteration', ...
        'ref_precision', 'rel_error', 'relerr_lt_u');
    fprintf(fid, '%8s  %22s  %9s  %10s  %13s  %14s  %12s\n', ...
        '------', '----------------', '-------', '---------', ...
        '-------------', '---------', '-----------');

    first_fail_exp = NaN;
    max_rel_error = -Inf;
    max_rel_error_exp = NaN;
    stopped_count = 0;
    inaccurate_count = 0;
    for exp_value = exponents
        rng(exp_value, 'twister');
        [p, ~, actual_condition] = GenPoly(degree, x, v, 10^exp_value);
        try
            [value, iter] = iter_CompHorner_accsum(p, x, max_iter);
        catch err
            if ~strcmp(err.identifier, 'iter_CompHorner_accsum:MaxIterations')
                rethrow(err);
            end
            if isnan(first_fail_exp)
                first_fail_exp = exp_value;
            end
            fprintf(fid, '1e%-6d  %22.15e  %9s  %10s  %13s  %14s  %12s\n', ...
                exp_value, actual_condition, 'false', 'fail', 'NaN', 'NaN', 'NaN');
            continue;
        end

        ref_precision = (iter + 1) * 53;
        [~, rel_error] = intlab_horner_relative_error(p, x, value, ref_precision);
        assert(isfinite(rel_error), 'Nonfinite relative error at degree %d, 1e%d.', degree, exp_value);
        accurate = rel_error < 2^-53;
        stopped_count = stopped_count + 1;
        inaccurate_count = inaccurate_count + ~accurate;
        if rel_error > max_rel_error
            max_rel_error = rel_error;
            max_rel_error_exp = exp_value;
        end
        accurate_text = 'false';
        if accurate
            accurate_text = 'true';
        end
        fprintf(fid, '1e%-6d  %22.15e  %9s  %10d  %13d  %14.6e  %12s\n', ...
            exp_value, actual_condition, 'true', iter, ref_precision, rel_error, accurate_text);
    end

    fprintf(fid, '\n');
    if isnan(first_fail_exp)
        fprintf(fid, 'First failure: none in 1e%d..1e%d.\n', exponents(1), exponents(end));
    else
        fprintf(fid, 'First failure: target condition 1e%d.\n', first_fail_exp);
    end
    fprintf(fid, 'Maximum relative error among stopped cases: %.6e at target condition 1e%d.\n', ...
        max_rel_error, max_rel_error_exp);
    fprintf('degree=%d: stopped=%d/%d, relerr>=u=%d, max_relerr=%.6e at 1e%d\n', ...
        degree, stopped_count, numel(exponents), inaccurate_count, max_rel_error, max_rel_error_exp);
    fprintf('Wrote: %s\n', output_file);
end
