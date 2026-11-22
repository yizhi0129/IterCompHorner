function [x, y] = two_product(a, b)
% error-free product: a * b = x + y
    x = a * b;

    [ah, al] = split_fp(a);
    [bh, bl] = split_fp(b);

    y = al * bl - (((x - ah * bh) - al * bh) - ah * bl);
end