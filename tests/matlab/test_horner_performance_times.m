function output_file = test_horner_performance_times()
%TEST_HORNER_PERFORMANCE_TIMES Compare direct times for two refined methods.

    root = fileparts(fileparts(mfilename('fullpath')));
    addpath(fullfile(root, 'src'));
    if exist('ufp_mex', 'file') ~= 3
        build_ufp_mex();
    end
    degree = 50;
    exponents = 5:5:30;
    trials = 10;
    seed = 0;
    max_iter = 100;
    min_abs = 2^-16;
    max_abs = 2^16;
    saved_rng = rng;
    cleanup_rng = onCleanup(@() rng(saved_rng)); %#ok<NASGU>
    rng(seed, 'twister');

    seconds = NaN(numel(exponents), trials, 2);
    actual_conditions = zeros(numel(exponents), trials);
    for c = 1:numel(exponents)
        for trial = 1:trials
            magnitudes = min_abs + (max_abs - min_abs) * rand(1, 2);
            inputs = (2 * randi([0, 1], 1, 2) - 1) .* magnitudes;
            x = inputs(1);
            v = inputs(2);
            [p, ~, actual_conditions(c, trial)] = GenPoly(degree, x, v, 10^exponents(c));

            methods = {@() iter_CompHorner_accsum(p, x, max_iter), ...
                       @() rpoly_eval_interval(p, x, max_iter, 80, "verified")};
            iter_ok = true;
            try
                methods{1}();
            catch err
                if ~strcmp(err.identifier, 'iter_CompHorner_accsum:MaxIterations')
                    rethrow(err);
                end
                iter_ok = false;
            end
            [~, ~, ~, rpoly_status] = rpoly_eval_interval(p, x, max_iter, 80, "verified");
            eligible = [iter_ok, strcmp(rpoly_status, 'No Error')];

            order = mod((0:1) + trial - 1, 2) + 1;
            for method = order(eligible(order))
                elapsed = timeit(methods{method}, 1);
                assert(isfinite(elapsed) && elapsed > 0, ...
                    'Invalid TIMEIT result for 1e%d, trial %d, method %d.', ...
                    exponents(c), trial, method);
                seconds(c, trial, method) = elapsed;
            end
        end
        fprintf('Finished target 1e%d (%d trials).\n', exponents(c), trials);
    end

    medians = reshape(median(seconds, 2, 'omitnan'), numel(exponents), 2);
    counts = reshape(sum(isfinite(seconds), 2), numel(exponents), 2);
    output_file = fullfile(root, 'results', 'horner_performance_times_degree50.txt');
    fid = fopen(output_file, 'w');
    assert(fid >= 0, 'Cannot open %s.', output_file);
    cleanup_file = onCleanup(@() fclose(fid)); %#ok<NASGU>
    fprintf(fid, 'Degree %d; %d independent trials per target; RNG seed %d (twister).\n', ...
        degree, trials, seed);
    fprintf(fid, 'x and v independently uniform on [-2^16,-2^-16] U [2^-16,2^16].\n');
    fprintf(fid, 'Methods: iter_CompHorner_accsum (iter_stop_test_accsum_mex), rpoly_eval_interval (verified, minimum 80 decimal digits).\n');
    fprintf(fid, 'Iteration cap: %d per method. Times are medians of per-trial TIMEIT estimates in seconds per evaluation, using identical p and x.\n', max_iter);
    fprintf(fid, 'Input generation and warmup excluded. Trials that reached the iteration cap without convergence are excluded per method.\n\n');
    fprintf(fid, 'target\tactual_cond_median\titer_CompHorner_accsum_seconds\trpoly_eval_interval_seconds\taccsum_n\trpoly_n\n');
    for c = 1:numel(exponents)
        fprintf(fid, '1e%d\t%.9e\t%.9e\t%.9e\t%d\t%d\n', exponents(c), ...
            median(actual_conditions(c, :)), medians(c, :), counts(c, :));
    end
    fprintf('Wrote: %s\n', output_file);
end
