function [tau, p_prime] = extract_vector(sigma, p)
% EFT extracting high order part of a vector
    
    tau = 0;
    n = length(p);
    q = zeros(1, n);
    p_prime = zeros(1, n);
    for i = 1:n
        q(i) = (sigma + p(i)) - sigma;
        p_prime(i) = p(i) - q(i);
        tau = tau + q(i);
    end    
end