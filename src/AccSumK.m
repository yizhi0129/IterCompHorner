function [Res, R, p_residual] = AccSumK(p, K)
%ACCSUMK Accurate summation computing K-fold faithful rounding.
%   Res = ACCSUMK(p, K) returns the K-term expansion produced by repeated
%   TransformK steps. Each step returns a faithfully rounded leading term
%   and advances the residual represented by R and p_residual.

    if ~isvector(p) || ~isreal(p)
        error('AccSumK:InvalidInput', 'p must be a real vector.');
    end
    if ~(isscalar(K) && isreal(K) && isfinite(K) && K >= 1 && K == floor(K))
        error('AccSumK:InvalidInput', 'K must be a positive integer.');
    end

    p_curr = p(:).';
    R_prev = 0;
    Res = zeros(1, K);

    for k = 1:K
        [Res(k), R_prev, p_curr] = TransformK(p_curr, R_prev);
        if is_in_U(Res(k))
            Res(k + 1:K) = 0;
            R = R_prev;
            p_residual = p_curr;
            return;
        end
    end

    R = R_prev;
    p_residual = p_curr;
end

function tf = is_in_U(x)
% U: subnormals, zero, and the two smallest normalized magnitudes.
    tf = isfinite(x) && abs(x) <= 2 * realmin;
end
