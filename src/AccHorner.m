function y = AccHorner(p, x, precision)
%ACCHORNER Evaluate a polynomial with MPFR Horner arithmetic.
%
%   y = AccHorner(p, x) evaluates p at x using 106-bit MPFR arithmetic
%   and rounds the final result back to double precision.
%
%   y = AccHorner(p, x, precision) uses the requested MPFR precision in
%   bits. The coefficient vector p follows MATLAB's decreasing-power order:
%   p = [a_n, a_{n-1}, ..., a_0].

    if nargin < 3
        precision = 106;
    end

    y = AccHorner_mex(p, x, precision);
end
