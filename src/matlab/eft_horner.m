function [res, p_pi, p_sigma] = eft_horner(p, x)
% error-free transformation of Horner scheme.
% p = [a_n, ..., a_0]
%
% Returns: (decreasing indices)
%   res      - Horner result
%   p_pi     - polynomial collecting multiplication errors
%   p_sigma  - polynomial collecting addition errors

    n = length(p) - 1;
    s = zeros(1, length(p));
    s(1) = p(1);
    p_pi = zeros(1, n);
    p_sigma = zeros(1, n);

    for k = 2:length(p)
        [prod_hi, prod_lo] = two_product(s(k-1), x);
        [s(k), sum_lo] = two_sum(prod_hi, p(k));

        p_pi(k-1) = prod_lo;
        p_sigma(k-1) = sum_lo;
    end

    res = s(length(p));
end