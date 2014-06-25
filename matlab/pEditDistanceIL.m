% pEditDistanceIL(a,b) treats a and b as internal loop sequences, ungapped, with * to separate the strands.  It returns several measures of edit distance.
% location == full
% d(1) = overall edit distance, minimum over both rotations
% d(2) = overall edit distance, original rotation
% d(3) = overall edit distance, rotated
% location == core / interior
% d(1) = core edit distance, minimum over both rotations
% d(2) = core edit distance, original order
% d(3) = core edit distance, rotated

function [d] = pEditDistanceIL(a,b,location)

d = zeros(1,3);                           % initial value

if nargin < 3,
  location = 'full';
end

i = strfind(a,'*');
j = strfind(b,'*');

if length(i) ~= 1 || length(j) ~= 1,
  d = Inf * ones(1,3);
elseif strcmpi(location,'full'),
  p = a(1:(i-1));                             % first sequence, first strand
  q = a((i+1):end);                           % first sequence, second strand
  r = b(1:(j-1));                             % second sequence, first strand
  s = b((j+1):end);                           % second sequence, second strand

  d(2) = EditDist(p,r)+EditDist(q,s);         % edit distance, original order
  d(3) = EditDist(p,s)+EditDist(q,r);         % edit distance, reversed order
  d(1) = min(d(2),d(3));                      % minimum over two orders
else
  p = a(2:(i-2));                             % first sequence, first strand
  q = a((i+2):(end-1));                       % first sequence, second strand
  r = b(2:(j-2));                             % second sequence, first strand
  s = b((j+2):(end-1));                       % second sequence, second strand

  d(2) = EditDist(p,r)+EditDist(q,s);         % edit distance, original order
  d(3) = EditDist(p,s)+EditDist(q,r);         % edit distance, reversed order
  d(1) = min(d(2),d(3));                      % minimum core distance
end

if 0 > 1,
  w = p([1 end]);                             % first, first, flank
  x = q([1 end]);                             % first, second, flank
  y = r([1 end]);                             % second, first, flank
  z = s([1 end]);                             % second, second, flank

  d(8) = EditDist(w,y)+EditDist(x,z);         % flank edit dist, original order
  d(9) = EditDist(w,z)+EditDist(x,y);         % flank edit dist, reversed order
  d(7) = min(d(5),d(6));                      % minimum flank distance
end

