x = Dual(2, 1);

f = func(x)
% exact derivative = 18.880984872278777

function [f] = func(x)
    f = exp(x)/sqrt(sin(x)^3 + cos(x)^3);
end
