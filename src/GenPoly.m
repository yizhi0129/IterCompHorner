function [p, vact, cact] = GenPoly(d, x, v, cexp, max_attempts)
% generate an ill-conditioned polynomial
% input:    d degree,
%           x argument,
%           v estimated p(x) value,
%           cexp expected conditioning number,
%           max_attempts optional retry limit;
% output:   p (Matlab decreasing order) coefficient vector,
%           vact p(x) evaluation,
%           cact conditioning number.

% add safety check: no underflow/overflow
    if nargin < 5
        max_attempts = 100;
    end

    for attempt = 1:max_attempts
        [p, vact, cact] = try_generate(d, x, v, cexp);

        if is_safe_polynomial(p, vact, cact)
            return;
        end
    end

    error('GenPoly:UnsafePolynomial', ...
          ['Failed to generate a polynomial without coefficient ', ...
           'underflow/overflow after %d attempts.'], max_attempts);
end

function [p, vact, cact] = try_generate(d, x, v, cexp)
    p = zeros(1, d + 1);
    d2 = ceil((d + 1) / 2);
    b = log2(cexp) + log2(abs(v));
    perm = randperm(d + 1);

    p(perm(1)) = (2 * rand - 1) / x^(d - perm(1) + 1);
    p(perm(2)) = (2 * rand - 1) * 2^b / x^(d - perm(2) + 1);
    for i = 3:d2
        p(perm(i)) = (2 * rand - 1) * 2^(b * rand) / ...
                     x^(d - perm(i) + 1);
    end

    log2v = log2(abs(v));
    m = (2 * rand(1, d - d2) - 1) .* 2.^linspace(b, log2v, d - d2);
    for i = d2 + 1:d
        p(perm(i)) = (m(i - d2) - AccHorner(p, x)) / ...
                     x^(d - perm(i) + 1);
    end

    p(perm(d + 1)) = (v - AccHorner(p, x)) / ...
                     x^(d - perm(d + 1) + 1);

    vact = AccHorner(p, x);
    cact = polyval(abs(p), abs(x)) / abs(vact);
end

function tf = is_safe_polynomial(p, vact, cact)
    tf = all(isfinite(p)) && ...
         all(p == 0 | abs(p) >= realmin) && ...
         isfinite(vact) && ...
         abs(vact) >= realmin && ...
         isfinite(cact);
end
