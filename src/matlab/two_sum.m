function [x, y] = two_sum(a, b)
% error-free sum: a + b = x + y 
    x = a + b;
    z = x - a;
    y = (a - (x - z)) + (b - z);
end