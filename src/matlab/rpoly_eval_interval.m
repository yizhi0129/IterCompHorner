function [value, value_int, iter, err, y] = rpoly_eval_interval(p, x, max_iter, long_precision, stop_mode)
%RPOLY_EVAL_INTERVAL Interval iterative refinement for Horner evaluation.
%   p is ordered as [a_n, a_{n-1}, ..., a_0].  Internally the C++ RPolyEval
%   indexing is used: coefficients and solution components are stored as
%   p_0,...,p_n and x_0,...,x_n, where x_0 = p(x).
%
%   value_int encloses p(x) for the converted double inputs. The loop
%   refines the Horner solution using enclosures of the residual correction.
%   long_precision is a minimum decimal working precision. The function
%   raises it as needed for exact local arithmetic.

    if nargin < 3
        max_iter = 1000;
    end
    if nargin < 4
        long_precision = 80;
    end
    if nargin < 5
        stop_mode = "verified";
    end

    validate_iter_inputs(p, max_iter, 'rpoly_eval_interval');
    stop_mode = validatestring(stop_mode, ["verified", "point", "either"], ...
                               'rpoly_eval_interval', 'stop_mode');

    if ~(isscalar(x) && isreal(x) && isfinite(x))
        error('rpoly_eval_interval:InvalidInput', ...
              'x must be a finite real scalar.');
    end

    if exist('intval', 'file') ~= 2
        error('rpoly_eval_interval:IntlabUnavailable', ...
              ['INTLAB is required. Run startintlab or add INTLAB to ', ...
               'the MATLAB path before calling this function.']);
    end
    if exist('long', 'file') ~= 2
        error('rpoly_eval_interval:IntlabUnavailable', ...
              'INTLAB long arithmetic is required for high-precision residuals.');
    end
    validate_long_precision(long_precision);

    p = full(double(p(:)));
    if any(~isfinite(p))
        error('rpoly_eval_interval:InvalidInput', 'p must have finite coefficients.');
    end
    p_asc = flipud(p);
    x = double(x);
    m = numel(p);
    n = m - 1;
    err = 'No Error';

    if n == 0
        value_int = intval(p(1));
        value = mid(value_int);
        iter = 0;
        y = p(1);
        return;
    end

    saved_precision = longprecision;
    cleanup_precision = onCleanup(@() longprecision(saved_precision)); %#ok<NASGU>
    % A product of two finite doubles lies on the 2^-2148 lattice and
    % has magnitude below 2^2048. Allow for all terms and guard digits.
    exact_bits = 4196 + ceil(log2(max_iter + 1)) + 34;
    longprecision(max(long_precision, ceil(exact_bits * log10(2))));
    p_long = long(p_asc);
    x_long = long(x);

    if n == 1
        exact_value = x_long * p_long(2) + p_long(1);
        value_int = infsup(round_down(exact_value), -round_down(-exact_value));
        value = mid(value_int);
        iter = 0;
        y = [p(1); value];
        return;
    end

    initial = initial_horner_solution_ascending(p_asc, x);
    if any(~isfinite(initial))
        initial = [zeros(n, 1); p_asc(end)];
    end
    approximation = long(initial);

    for k = 1:max_iter
        [correction, value_int] = solve_correction(p_long, x_long, x, approximation);
        correction_mid = mid(infsup(correction(:, 1), correction(:, 2)));
        current_value = double(approximation(1));
        value = mid(value_int);
        lower = inf(value_int);
        upper = sup(value_int);
        finite_bounds = isfinite(lower) && isfinite(upper);
        point_stop = finite_bounds && isfinite(current_value) && ...
            abs(correction_mid(1)) <= eps(abs(current_value));
        verified_stop = finite_bounds && ...
            long(upper) - long(lower) <= long(eps(value));

        approximation = approximation + long(correction_mid);
        stop_success = should_stop(verified_stop, point_stop, stop_mode);
        if stop_success || k == max_iter
            iter = k;
            y = flipud(double(approximation));
            y(end) = value;
            break;
        end
    end

    if ~stop_success
        err = 'Iteration failed';
    end
end

function tf = should_stop(verified_stop, point_stop, stop_mode)
    switch string(stop_mode)
        case "verified"
            tf = verified_stop;
        case "point"
            tf = point_stop;
        case "either"
            tf = verified_stop || point_stop;
    end
end

function validate_long_precision(long_precision)
    if ~(isscalar(long_precision) && isreal(long_precision) && ...
         isfinite(long_precision) && long_precision >= 1 && ...
         long_precision == floor(long_precision))
        error('rpoly_eval_interval:InvalidInput', ...
              'long_precision must be a positive integer.');
    end
end

function y = initial_horner_solution_ascending(p, x)
    m = numel(p);
    n = m - 1;
    y = zeros(m, 1);
    y(m) = p(m);

    for idx = n:-1:1
        y(idx) = y(idx + 1) * x + p(idx);
    end
end

function [correction, value_int] = solve_correction(p, x_long, x, approximation)
    m = size(approximation, 1);
    correction = zeros(m, 2);
    for idx = m - 1:-1:1
        next = correction(idx + 1, :);
        if x < 0
            next = next([2, 1]);
        end
        if any(~isfinite(next))
            correction(1:idx, 1) = -Inf;
            correction(1:idx, 2) = Inf;
            value_int = infsup(-Inf, Inf);
            return;
        end

        residual = p(idx) - approximation(idx) + x_long * approximation(idx + 1);
        lower = residual + x_long * long(next(1));
        upper = residual + x_long * long(next(2));
        correction(idx, :) = [round_down(lower), -round_down(-upper)];
    end
    value_int = infsup(round_down(approximation(1) + lower), ...
                       -round_down(-approximation(1) - upper));
end

function bound = round_down(exact)
    % Use the conversion only as a starting point. Compare in exact long
    % arithmetic before accepting the endpoint, independent of setround.
    bound = min(double(exact), realmax);
    while isfinite(bound) && long(bound) > exact
        if bound == 0
            bound = -eps(0);
        else
            bits = typecast(bound, 'uint64');
            if bound < 0
                bits = bits + uint64(1);
            else
                bits = bits - uint64(1);
            end
            bound = typecast(bits, 'double');
        end
    end
end
