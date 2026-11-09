function y = ufp(x)
%UFP Unit in the first place of a floating-point value.
%   y = UFP(x) returns 2^floor(log2(abs(x))) elementwise.

    y = 2 .^ floor(log2(abs(x)));
    y(x == 0) = 0;
end
