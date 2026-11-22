function [tf, lhs, rhs] = iter_stop_test_accsum_mex(s, k, theta, E, gamma)
%ITER_STOP_TEST_ACCSUM_MEX Revised Algorithm 13 with scalar error bound.
%   theta >= 0 bounds abs(exact_sum(H)-s); it is NOT a theta_list.
%   k is retained for the Algorithm 13 interface; theta is already bounded,
%   so no additional division by 1-(k-1)*u is made here.
%   A zero or underflowed RHS cannot satisfy the strict inequality.
    if ~(isscalar(s) && isreal(s) && isfinite(s) && ...
         isscalar(k) && isfinite(k) && k >= 1 && k == floor(k) && ...
         isscalar(theta) && isreal(theta) && isfinite(theta) && theta >= 0 && ...
         isscalar(E) && isreal(E) && isfinite(E) && E >= 0 && ...
         isscalar(gamma) && isreal(gamma) && isfinite(gamma) && gamma >= 0)
        error('iter_stop_test_accsum_mex:InvalidInput', ...
            'Require finite real scalar inputs, k>=1 integer, and theta,E,gamma>=0.');
    end
    u = 2^-53;
    b = (1+2*u)*E*(1+gamma/(2*(1-2*u)))/(1-4*u);
    lhs = (theta+b)/(1-2*u);
    rhs = u*ufp_mex(s);
    tf = isfinite(lhs) && lhs < rhs;
end
