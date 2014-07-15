% pEditDistanceIL(a,b) treats a and b as a hairpin loop sequence.  It returns full or interior/core edit distance, as requested.

function [d] = pEditDistanceHL(a,b,location)

if nargin < 3,
	location = 'full';
end

d = 0;                           % initial value

i = strfind(a,'*');
j = strfind(b,'*');

if length(i) ~= 0 || length(j) ~= 0,
  d = Inf;
elseif strcmpi(location,'full'),
  d = EditDist(a,b);
else
  w = a(2:(end-1));                           % first, core
  x = b(2:(end-1));                           % second, core
  d = EditDist(w,x);
end
