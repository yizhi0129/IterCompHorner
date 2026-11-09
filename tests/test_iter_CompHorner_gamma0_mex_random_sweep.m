function output_file = test_iter_CompHorner_gamma0_mex_random_sweep()
%TEST_ITER_COMPHORNER_GAMMA0_MEX_RANDOM_SWEEP Test random x and v at degree 50.

    if exist('ufp_mex', 'file') ~= 3
        build_ufp_mex();
    end

    degree = 50;
    exponents = 2:35;
    repeats = 20;
    max_iter = 100;
    seed = 0;
    min_abs = 2^-16;
    max_abs = 2^16;
    saved_rng = rng;
    cleanup_rng = onCleanup(@() rng(saved_rng)); %#ok<NASGU>
    rng(seed, 'twister');

    root = fileparts(mfilename('fullpath'));
    output_file = fullfile(root, ...
        'iter_CompHorner_gamma0_mex_degree50_random_x_v_sweep_1e2_1e35_relerr_output.txt');
    fid = fopen(output_file, 'w');
    assert(fid >= 0, 'Cannot open %s.', output_file);
    cleanup_file = onCleanup(@() fclose(fid)); %#ok<NASGU>

    fprintf(fid, '=== iter_CompHorner AccSumK-H / MEX-stop sweep with random x and v ===\n\n');
    fprintf(fid, 'implementation = iter_CompHorner_accsum\nstop = iter_stop_test_accsum_mex\n');
    fprintf(fid, 'degree = %d\nmax_iter = %d\nrepeats_per_target = %d\nseed = %d (twister)\n', ...
        degree, max_iter, repeats, seed);
    fprintf(fid, 'x and v are independent uniform draws on [-2^16, -2^-16] U [2^-16, 2^16].\n\n');
    fprintf(fid, '%24s  %24s  %8s  %22s  %9s  %10s  %13s  %22s  %12s\n', ...
        'x', 'v', 'target', 'actual_condition', 'stopped', 'iteration', ...
        'ref_precision', 'rel_error', 'relerr_lt_u');

    stopped_count = 0;
    inaccurate_count = 0;
    max_rel_error = NaN;
    for exponent = exponents
        for trial = 1:repeats
            magnitudes = min_abs + (max_abs - min_abs) * rand(1, 2);
            signs = 2 * randi([0, 1], 1, 2) - 1;
            x = signs(1) * magnitudes(1);
            v = signs(2) * magnitudes(2);
            [p, ~, actual_condition] = GenPoly(degree, x, v, 10^exponent);

            try
                [value, iter] = iter_CompHorner_accsum(p, x, max_iter);
            catch err
                if ~strcmp(err.identifier, 'iter_CompHorner_accsum:MaxIterations')
                    rethrow(err);
                end
                fprintf(fid, '%24.17g  %24.17g  1e%-6d  %22.15e  %9s  %10s  %13s  %22s  %12s\n', ...
                    x, v, exponent, actual_condition, 'false', 'fail', 'NaN', 'NaN', 'NaN');
                continue;
            end

            ref_precision = (iter + 1) * 53;
            [~, rel_error] = intlab_horner_relative_error(p, x, value, ref_precision);
            assert(isfinite(rel_error), 'Nonfinite relative error at 1e%d, trial %d.', exponent, trial);
            accurate = rel_error < 2^-53;
            stopped_count = stopped_count + 1;
            inaccurate_count = inaccurate_count + ~accurate;
            max_rel_error = max(max_rel_error, rel_error);
            accuracy_labels = {'false', 'true'};
            fprintf(fid, '%24.17g  %24.17g  1e%-6d  %22.15e  %9s  %10d  %13d  %22.16e  %12s\n', ...
                x, v, exponent, actual_condition, 'true', iter, ref_precision, ...
                rel_error, accuracy_labels{1 + accurate});
        end
    end

    fprintf(fid, '\nStopped: %d/%d\nStopped with relerr >= u: %d\n', ...
        stopped_count, numel(exponents) * repeats, inaccurate_count);
    fprintf(fid, 'Maximum relative error among stopped cases: %.16e\n', max_rel_error);
end
