function [x, y] = vect_two_sum(a, b)
% error-free elementwise sum: a + b = x + y
% different length: align first elements and pad zeros to the end.
    if ~isvector(a) || ~isvector(b)
        error('vect_two_sum:InvalidInput', 'Inputs must be vectors.');
    end

    n = max(numel(a), numel(b));
    use_column = (iscolumn(a) && numel(a) >= numel(b)) || ...
                 (iscolumn(b) && numel(b) > numel(a));

    if use_column
        aa = zeros(n, 1);
        bb = zeros(n, 1);
    else
        aa = zeros(1, n);
        bb = zeros(1, n);
    end

    aa(1:numel(a)) = a;
    bb(1:numel(b)) = b;

    x = aa + bb;
    z = x - aa;
    y = (aa - (x - z)) + (bb - z);
end
