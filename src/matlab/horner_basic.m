function y = horner_basic(p, x)
% p = [a_n, a_{n-1}, ..., a_0]
    y = p(1);
    for k = 2:length(p)
        y = y * x + p(k);
    end
end