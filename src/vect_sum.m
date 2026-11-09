function x = vect_sum(a, b)
% add two vectors aligned at the first element, padding zeros at end.
    if ~isvector(a) || ~isvector(b)
        error('vect_sum:InvalidInput', 'Inputs must be vectors.');
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
end
