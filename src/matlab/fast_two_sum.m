function [x, y] = fast_two_sum(a, b)
% a + b = x + y with |a| ≥ |b|
    
    x = a + b;
    q = x - a;
    y = b - q;
end