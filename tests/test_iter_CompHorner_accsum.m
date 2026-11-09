function summary = test_iter_CompHorner_accsum()
%TEST_ITER_COMPHORNER_ACCSUM Check revised algorithm against INTLAB.
%   Test the 156 original sweep cases, seed-30 known failure, and all
%   failed cases saved by the latest failure-frequency experiment.
    root = fileparts(mfilename('fullpath'));
    if exist('ufp_mex', 'file') ~= 3
        build_ufp_mex();
    end
    saved_rng = rng;
    cleanup_rng = onCleanup(@() rng(saved_rng)); %#ok<NASGU>
    degrees = [10,40,70,100];
    cases = zeros(0, 3);
    for d = degrees
        cases = [cases; [repmat(d,39,1), (2:40).', (2:40).']]; %#ok<AGROW>
    end
    data_file = fullfile(root, 'gamma0_mex_failure_frequency_results.mat');
    if exist(data_file, 'file') == 2
        data = load(data_file, 'results');
        cases = [cases; data.results.failed_cases(:,1:3)];
    end
    cases = unique(cases, 'rows', 'stable');

    assert(iter_CompHorner_accsum(3, 2, 1) == 3);
    % The strict stop inequality cannot certify cancellation to zero.
    expect_error(@() iter_CompHorner_accsum([1,-1],1,1), ...
        'iter_CompHorner_accsum:MaxIterations');
    assert(iter_CompHorner_accsum([2,3], 0, 1) == 3);
    [~, k] = iter_CompHorner_accsum([1,1], 1, 1);
    assert(k == 1);
    assert(~iter_stop_test_accsum_mex(0,1,0,0,0));
    assert(~iter_stop_test_accsum_mex(realmin,1,0,0,0));
    expect_error(@() iter_stop_test_accsum_mex(1,2,-1,0,0), ...
        'iter_stop_test_accsum_mex:InvalidInput');
    rng(30,'twister');
    difficult = GenPoly(100,1.03,1e-3,1e30);
    expect_error(@() iter_CompHorner_accsum(difficult,1.03,1), ...
        'iter_CompHorner_accsum:MaxIterations');

    output = fullfile(root,'iter_CompHorner_accsum_validation_output.txt');
    fid = fopen(output,'w');
    assert(fid >= 0);
    cleanup_file = onCleanup(@() fclose(fid)); %#ok<NASGU>
    fprintf(fid,'Revised iter_CompHorner: AccSumK(H,k), scalar-theta MEX stop\n');
    fprintf(fid,'x=1.03 v=1e-3 max_iter=100; INTLAB reference=512 bits\n');
    fprintf(fid,'degree target_exp seed actual_condition stopped iteration rel_error relerr_lt_u\n');
    stopped_count = 0;
    inaccurate = 0;
    max_error = 0;
    for j = 1:size(cases,1)
        d = cases(j,1); e = cases(j,2); seed = cases(j,3);
        rng(seed,'twister');
        [p,~,actual] = GenPoly(d,1.03,1e-3,10^e);
        try
            [s,k,trace] = iter_CompHorner_accsum(p,1.03,100);
        catch err
            if ~strcmp(err.identifier,'iter_CompHorner_accsum:MaxIterations')
                rethrow(err);
            end
            fprintf(fid,'%d %d %d %.17e false fail NaN NaN\n',d,e,seed,actual);
            continue;
        end
        assert(size(trace,1) == k && all(trace(:,4) >= 0));
        [~,relerr] = intlab_horner_relative_error(p,1.03,s,512);
        % Independently verify the scalar bound on each accumulated H sum.
        exact_sum = long(0);
        for t = 1:k
            exact_sum = exact_sum + long(trace(t,3));
            assert(abs(exact_sum-long(trace(t,2))) <= long(trace(t,4)), ...
                'theta failed to bound sum(H)-s at degree=%d, exponent=%d, seed=%d, k=%d.', ...
                d,e,seed,t);
        end
        accurate = isfinite(relerr) && relerr < 2^-53;
        stopped_count = stopped_count+1;
        inaccurate = inaccurate+~accurate;
        max_error = max(max_error,relerr);
        fprintf(fid,'%d %d %d %.17e true %d %.17e %d\n',d,e,seed,actual,k,relerr,accurate);
    end
    summary = struct('tested',size(cases,1),'stopped',stopped_count, ...
        'failed',size(cases,1)-stopped_count,'inaccurate',inaccurate, ...
        'max_rel_error',max_error);
    for dest = [1,fid]
        fprintf(dest,'Tested=%d stopped=%d failed=%d relerr>=u=%d max_relerr=%.17e\n', ...
            summary.tested,summary.stopped,summary.failed,inaccurate,max_error);
        fprintf(dest,'All recorded theta bounds passed high-precision checks of abs(sum(H)-s).\n');
    end
end

function expect_error(fun,id)
    try
        fun();
    catch err
        assert(strcmp(err.identifier,id),'Unexpected error: %s',err.identifier);
        return;
    end
    error('test_iter_CompHorner_accsum:ExpectedError','Expected %s.',id);
end
