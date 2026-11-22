function [xh, xl] = split_fp(a)
% split double
    c = (2^27 + 1) * a;
    xh = c - (c - a);
    xl = a - xh;
end