%TEST_ITER_COMPHORNER_CONDITION_SWEEP_1E5_1E40_UNDERFLOW
%  Sweep target condition numbers 1e5, 1e6, ..., 1e40.  For each generated
%  polynomial, run the iter_CompHorner recurrence, record the stopping
%  iteration, and count subnormal nonzero intermediate values.

clc;

degree = 10;
x = 1.03;
v = 1e-3;
max_iter = 1000;
exponents = 5:40;

fprintf('=== iter_CompHorner condition sweep with underflow diagnostics ===\n\n');
fprintf('degree = %d\n', degree);
fprintf('x = %.17g\n', x);
fprintf('v = %.17g\n', v);
fprintf('max_iter = %d\n', max_iter);
fprintf('underflow count = number of nonzero intermediate values with abs(value) < realmin\n\n');

fprintf('%8s  %22s  %9s  %10s  %14s  %12s  %14s\n', ...
        'target', 'actual_condition', 'stopped', 'iteration', ...
        'rel_error', 'underflows', 'min_nonzero_abs');
fprintf('%8s  %22s  %9s  %10s  %14s  %12s  %14s\n', ...
        '------', '----------------', '-------', '---------', ...
        '---------', '----------', '---------------');

first_fail_exp = NaN;
total_underflows = 0;

for exp_value = exponents
    rng(exp_value);
    target_condition = 10^exp_value;
    [p, ref, actual_condition] = GenPoly(degree, x, v, target_condition);

    result = run_iter_with_underflow_count(p, x, max_iter);
    total_underflows = total_underflows + result.underflow_count;

    if result.stopped
        rel_error = abs(result.value - ref) / max(abs(ref), realmin);
        fprintf('1e%-6d  %22.15e  %9s  %10d  %14.6e  %12d  %14.6e\n', ...
                exp_value, actual_condition, 'true', result.iter, ...
                rel_error, result.underflow_count, result.min_nonzero_abs);
    else
        if isnan(first_fail_exp)
            first_fail_exp = exp_value;
        end
        fprintf('1e%-6d  %22.15e  %9s  %10s  %14s  %12d  %14.6e\n', ...
                exp_value, actual_condition, 'false', 'fail', 'NaN', ...
                result.underflow_count, result.min_nonzero_abs);
    end
end

fprintf('\n');
if isnan(first_fail_exp)
    fprintf('First failure: none in 1e5..1e40.\n');
else
    fprintf('First failure: target condition 1e%d.\n', first_fail_exp);
end
fprintf('Total underflow count across sweep: %d.\n', total_underflows);

function result = run_iter_with_underflow_count(p, x, max_iter)
    validate_iter_inputs(p, max_iter, 'iter_CompHorner');

    p = p(:).';
    d = max(numel(p) - 2, 0);
    u = 2^-53;
    gamma = (4 * d * u) / (1 - 4 * d * u);

    underflow_count = 0;
    min_nonzero_abs = Inf;

    [h, p_pi, p_sigma] = eft_horner(p, x);
    [underflow_count, min_nonzero_abs] = update_underflow_stats( ...
        underflow_count, min_nonzero_abs, h, p_pi, p_sigma);

    e_family = accsumk_vector_full({p_pi, p_sigma});
    [underflow_count, min_nonzero_abs] = update_underflow_stats( ...
        underflow_count, min_nonzero_abs, e_family{:});

    res = h;
    E = horner_basic(abs(e_family{1}), abs(x));
    [underflow_count, min_nonzero_abs] = update_underflow_stats( ...
        underflow_count, min_nonzero_abs, res, E);

    if iter_stop_test_gamma(res, 1, [], E, gamma)
        result = make_result(true, 1, res, underflow_count, min_nonzero_abs);
        return;
    end

    theta_list = zeros(1, max_iter - 1);
    for k = 2:max_iter
        [h, p_pi, p_sigma] = eft_horner(e_family{1}, x);
        [underflow_count, min_nonzero_abs] = update_underflow_stats( ...
            underflow_count, min_nonzero_abs, h, p_pi, p_sigma);

        e_family = accsumk_vector_full([{p_pi, p_sigma}, e_family(2:end)]);
        [underflow_count, min_nonzero_abs] = update_underflow_stats( ...
            underflow_count, min_nonzero_abs, e_family{:});

        [res, theta_k] = two_sum(res, h);
        theta_list(k - 1) = theta_k;

        E = horner_basic(abs(e_family{1}), abs(x));
        [underflow_count, min_nonzero_abs] = update_underflow_stats( ...
            underflow_count, min_nonzero_abs, res, theta_k, E);

        if iter_stop_test_gamma(res, k, theta_list(1:k - 1), E, gamma)
            result = make_result(true, k, res, underflow_count, min_nonzero_abs);
            return;
        end
    end

    result = make_result(false, max_iter, res, underflow_count, min_nonzero_abs);
end

function result = make_result(stopped, iter, value, underflow_count, min_nonzero_abs)
    if isinf(min_nonzero_abs)
        min_nonzero_abs = NaN;
    end

    result = struct( ...
        'stopped', stopped, ...
        'iter', iter, ...
        'value', value, ...
        'underflow_count', underflow_count, ...
        'min_nonzero_abs', min_nonzero_abs);
end

function [underflow_count, min_nonzero_abs] = update_underflow_stats( ...
        underflow_count, min_nonzero_abs, varargin)
    for idx = 1:numel(varargin)
        values = varargin{idx};
        values = values(:);
        finite_nonzero = values(isfinite(values) & values ~= 0);
        if ~isempty(finite_nonzero)
            min_nonzero_abs = min(min_nonzero_abs, min(abs(finite_nonzero)));
            underflow_count = underflow_count + ...
                              nnz(abs(finite_nonzero) < realmin);
        end
    end
end
