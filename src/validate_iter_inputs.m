function validate_iter_inputs(p, max_iter, error_prefix)
    if ~isvector(p) || ~isreal(p)
        error([error_prefix ':InvalidInput'], 'p must be a real vector.');
    end

    if isempty(p)
        error([error_prefix ':InvalidInput'], 'p must not be empty.');
    end

    if ~(isscalar(max_iter) && isreal(max_iter) && isfinite(max_iter) && ...
         max_iter >= 1 && max_iter == floor(max_iter))
        error([error_prefix ':InvalidInput'], ...
              'max_iter must be a positive integer.');
    end
end
