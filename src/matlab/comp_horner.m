function [res, h, c, p_pi, p_sigma] = comp_horner(p, x)
% Compensated Horner evaluation
% p = [a_n, ..., a_0]

    [h, p_pi, p_sigma] = eft_horner(p, x);

    p_r = p_pi + p_sigma;
    c = horner_basic(p_r, x);
    res = h + c;
end
