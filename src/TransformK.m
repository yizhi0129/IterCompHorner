function [res, R, p_prime] = TransformK(p, rho)
    u = 2^-53;
    if ~isvector(p) || ~isreal(p)
        error('TransformK:InvalidInput', 'p must be a real vector.');
    end
    if ~(isscalar(rho) && isreal(rho) && isfinite(rho))
        error('TransformK:InvalidInput', 'rho must be a finite real scalar.');
    end

    M = ceil(log2(length(p) + 2));
    Phi = 2^(2 * M) * u;

    [tau1, tau2, p_prime] = Transform(p, rho, Phi);

    res = tau1 + (tau2 + sum(p_prime));
    R = tau2 - (res - tau1);
end
