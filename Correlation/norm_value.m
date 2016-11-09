function [ p ] = norm_value(test_value,base_value, sigma)
%
    p = exp(-0.5 * ((test_value - base_value)/sigma)^2)/(sqrt(2*pi) * sigma);
end

