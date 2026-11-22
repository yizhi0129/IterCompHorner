function [tau1, tau2, p_m, sigma, M] = Transform(p0, rho, Phi)
% vector p0 plus offset rho
% default Phi = 2
% positive parameter Phi
    u = 2^-53;
    if nargin < 3
        Phi = 2;
    end

    if ~isvector(p0) || ~isreal(p0)
        error('Transform:InvalidInput', 'p0 must be a real vector.');
    end
    if ~(isscalar(rho) && isreal(rho) && isfinite(rho))
        error('Transform:InvalidInput', 'rho must be a finite real scalar.');
    end
    if ~(isscalar(Phi) && isreal(Phi) && isfinite(Phi) && Phi > 0)
        error('Transform:InvalidInput', 'Phi must be a positive finite scalar.');
    end

    mu = max(abs(p0));
    if mu == 0
        tau1 = rho;
        tau2 = 0;
        p_m = zeros(size(p0));
        sigma = 0;
        M = 0;
        return;
    end

    M = ceil(log2(length(p0) + 2));
    sigma_prev = 2^(M + ceil(log2(mu)));

    t_prev = rho;
    p_prev = p0;
    eta = realmin * u;
    sigma_min = 0.5 * u^-1 * eta;

    while true
        [tau_m, p_m] = extract_vector(sigma_prev, p_prev);
        t_m = t_prev + tau_m;
        sigma_m = (2^M * u) * sigma_prev;

        if abs(t_m) >= Phi * sigma_prev || sigma_prev <= sigma_min
            sigma = sigma_prev;
            [tau1, tau2] = two_sum(t_prev, tau_m);
            return;
        end

        t_prev = t_m;
        p_prev = p_m;
        sigma_prev = sigma_m;
    end
end
