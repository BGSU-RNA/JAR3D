% pEditDistance calculates the pairwise edit distances between the sequences in F and G, which are assumed to be of type loopType.  It returns a length(F) x length(G) matrix of pairwise edit distances.
% For internal loops, D1 is the minimum edit distance over both strand orders, D2 is the edit distance in the given order, D3 reversed order
% location can be 'full' or 'core'

function [D1,D2,D3] = pEditDistance(F,G,loopType,location)

m = length(F);

if isempty(G),
  G = F;
  self = 1;
else
  self = 0;
end

n = length(G);

if nargin < 3,
  loopType = 'HL';
end

if nargin < 4,
  location = 'Full';
end

D1 = Inf * ones(m,n);
D2 = Inf * ones(m,n);
D3 = Inf * ones(m,n);

if strcmp(class(F),'struct') && strcmp(class(G),'struct'),
  F = lower(F);
  G = lower(G);
  if strcmp(loopType,'HL'),
    for a = 1:m,
      if self > 0,
        D2(a,a) = 0;
        D3(a,a) = 0;
        for b = (a+1):n,
          d = pEditDistanceHL(F(a).Sequence,G(b).Sequence,location);
          D2(a,b) = d(1);
          D3(a,b) = d(1);
          D2(b,a) = d(1);
          D3(b,a) = d(1);
        end
      else
        for b = 1:n,
          d = pEditDistanceHL(F(a).Sequence,G(b).Sequence,location);
          D1(a,b) = d(1);
          D2(a,b) = d(1);
          D3(a,b) = d(1);
        end
      end
    end
  end
  if strcmp(loopType,'IL'),
    for a = 1:m,
      if self > 0,
        for b = a:n,
          d = pEditDistanceIL(F(a).Sequence,G(b).Sequence,location);
          D2(a,b) = d(2);
          D3(a,b) = d(3);
          D2(b,a) = d(2);
          D3(b,a) = d(3);
        end
      else
        for b = 1:n,
          d = pEditDistanceIL(F(a).Sequence,G(b).Sequence,location);
          D2(a,b) = d(2);
          D3(a,b) = d(3);
        end
      end
    end
    if sum(sum(D2)) < sum(sum(D3)),
      D1 = D2;
    else
      D1 = D3;
    end
  end
end

