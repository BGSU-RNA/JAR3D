% use the probabilities in x to generate a matrix with M rows and N columns

function [r] = rando(x,M,N)

n = length(x);

if nargin == 2,
	N = 1;
end

if nargin == 1,
	M = 1;
	N = 1;
end

r = zeros(M,N);

for a = 1:M,
	for b = 1:N,
		u = rand;
		i = 1;
		s = x(1);

		while ((u > s) & (i < n)),
		  i=i+1;
		  s=s+x(i);
		end

		r(a,b)=i;
	end
end
