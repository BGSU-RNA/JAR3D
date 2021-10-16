% ---------------------------------- Poisson distribution for lengths -------

function [d] = subPoisson(m)

n = max(3,2*m);

d = exp(-m) * (m .^ (0:n)) ./ factorial(0:n);

d = d / sum(d);                     % normalize
