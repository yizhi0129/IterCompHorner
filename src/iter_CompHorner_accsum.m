function [res, i, trace] = iter_CompHorner_accsum(p, x, max_iter)
%ITER_COMPHORNER_ACCSUM Iterative Horner with accurate resummation of H.
%   [res,i] = iter_CompHorner_accsum(p,x,max_iter) implements the revised
%   Algorithm 12. At k=2 use TwoSum; at k>=3 use AccSumK(H(1:k),k).
%   The stop test takes a nonnegative scalar summation-error bound, rather
%   than the old list of errors from successive additions.
%   Optional trace records [k,s,h,theta,E,lhs,rhs,h2] for each iteration.
%
%   Error bounds inherit the EFT/AccSumK assumptions: round-to-nearest,
%   no overflow, and validity of their residual identities. In particular,
%   the h2 bound assumes h2 faithfully rounds the exact sum(H)-res.
%   This implementation does not establish a new underflow error theorem.

    if nargin < 3
        max_iter = 10;
    end
    validate_iter_inputs(p, max_iter, 'iter_CompHorner_accsum');
    if ~isa(p, 'double') || issparse(p) || any(~isfinite(p)) || ...
            ~isa(x, 'double') || ~isscalar(x) || ~isreal(x) || ~isfinite(x)
        error('iter_CompHorner_accsum:InvalidInput', ...
            'p must be a full finite real double vector; x a finite real double scalar.');
    end
    p = p(:).';
    trace = zeros(0, 8);
    % These cases can be evaluated exactly without any error transforms.
    if numel(p) == 1 || x == 0 || all(p == 0)
        res = p(end);
        i = 1;
        return;
    end
    u = 2^-53;
    M = ceil(log2(max_iter+2));
    if 2^(2*M)*u > 1
        error('iter_CompHorner_accsum:MaxIterationsTooLarge', ...
            'max_iter violates the AccSumK precondition 2^(2*M)*u <= 1.');
    end
    d = max(numel(p)-2, 0);
    gamma = (4*d*u)/(1-4*d*u);
    H = zeros(1, max_iter);

    for k = 1:max_iter
        if k == 1
            [h, p_pi, p_sigma] = eft_horner(p, x);
            e_family = accsumk_vector_full({p_pi, p_sigma});
        else
            [h, p_pi, p_sigma] = eft_horner(e_family{1}, x);
            e_family = accsumk_vector_full([{p_pi, p_sigma}, e_family(2:end)]);
        end
        H(k) = h;
        h2 = 0;
        if k == 1
            res = h;
            theta = 0;
        elseif k == 2
            [res, residual] = two_sum(H(1), H(2));
            theta = abs(residual);
            h2 = residual;
        else
            expansion = AccSumK(H(1:k), k);
            res = expansion(1);
            h2 = expansion(2);
            theta = abs(h2)*(1+2*u)/(1-2*u);
        end
        E = horner_basic(abs(e_family{1}), abs(x));
        if ~isfinite(res) || ~isfinite(theta) || ~isfinite(E)
            error('iter_CompHorner_accsum:NonfiniteIntermediate', ...
                'Nonfinite intermediate at iteration %d.', k);
        end
        [stopped, lhs, rhs] = iter_stop_test_accsum_mex(res, k, theta, E, gamma);
        if nargout >= 3
            trace(k, :) = [k, res, h, theta, E, lhs, rhs, h2]; %#ok<AGROW>
        end
        if stopped
            i = k;
            return;
        end
    end
    error('iter_CompHorner_accsum:MaxIterations', ...
        'Stopping criterion was not met within %d iterations.', max_iter);
end
