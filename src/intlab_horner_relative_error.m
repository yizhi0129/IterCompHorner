function [ref_value, rel_error] = intlab_horner_relative_error(p, x, value, precision)
%INTLAB_HORNER_RELATIVE_ERROR  High-precision INTLAB reference and rel error.
%
%   [ref_value, rel_error] = intlab_horner_relative_error(p, x, value, precision)
%   evaluates p(x) with INTLAB long arithmetic at the requested precision,
%   then computes abs(value - ref) / max(abs(ref), realmin) at that same
%   precision.  p is in decreasing-power order.

    if exist('long', 'file') ~= 2
        if exist('startintlab', 'file') == 2
            startintlab;
        end
    end

    if exist('long', 'file') ~= 2
        error('intlab_horner_relative_error:INTLABUnavailable', ...
              'INTLAB long arithmetic is required.');
    end

    if ~(isscalar(precision) && isreal(precision) && isfinite(precision) && ...
         precision >= 1 && precision == floor(precision))
        error('intlab_horner_relative_error:InvalidPrecision', ...
              'precision must be a positive integer.');
    end

    longprecision(precision);

    p = p(:).';
    ncoeff = numel(p);
    p_long = long(p);
    x_long = long(x);
    value_long = long(value);

    ref_long = p_long(1);
    for idx = 2:ncoeff
        ref_long = ref_long * x_long + p_long(idx);
    end

    denom = abs(ref_long);
    if denom < long(realmin)
        denom = long(realmin);
    end

    rel_long = abs(value_long - ref_long) / denom;

    ref_value = double(ref_long);
    rel_error = double(rel_long);
end
