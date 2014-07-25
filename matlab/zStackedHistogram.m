% zStackedHistogram histograms the data points in X but colors the bars according to the values in Y.  Y needs to have integer values

function [counts,bins] = zStackedHistogram(X,Y,percentage,nbins,orientation)

if nargin < 2,
  Y = ones(size(X));
end

if length(Y) == 1,
  Y = Y * ones(size(X));
end

if length(Y) ~= length(X),
  fprintf('zStackedHistogram objects that X and Y have different lengths\n');
end

if nargin < 3,
  percentage = 0;
end

if nargin < 4,
  nbins = 30;
  if all(X == round(X)),               % integer values
    nbins = range(X);
  end
end

if nargin < 5,
  orientation = 'h';
end

if length(X) > 0,

  if 0 > 1,
    score = log(X);                    % put scores on a log scale
  else
    score = X;
  end

  if length(nbins) > 1,
    B = nbins;
  else
    B = min(score):range(score)/nbins:max(score);
  end

  if all(score == round(score)) && range(score) == nbins,% histograming integers
    B = (min(score)-0.5):(max(score)+0.5);
  end

  clear counts
  clear bins

  for d = min(Y):max(Y),
    i = find(Y == d);
    nc = histc(score(i),B);  % bins is same as B
    nc(end-1) = nc(end-1) + nc(end);
    nc = nc(1:(end-1));
    counts(:,d-min(Y)+1) = nc;  % bins is same as B
  end
  
  [s,t] = size(counts);

  counts = [zeros(s,1) counts];

  counts = cumsum(counts')';

  bins = B;

  if percentage > 0,
    zStackedBar(B,counts/length(X),[],orientation)
  else
    zStackedBar(B,counts,[],orientation)
  end

else
  counts = [];
  bins = [];
end

%  counts = 2.2 * counts / max(max(counts));
%  counts = 0.1*exp(counts);
%  counts = counts/max(max(counts));

% clf; [counts,bins] = zStackedHistogram(rand(1000,1),ceil(3*rand(1000,1)),1,10); sum(sum(counts(:,end)))

% clf; [counts,bins] = zStackedHistogram(rand(1000,1),ceil(3*rand(1000,1)),1,0:0.1:10); sum(sum(counts(:,end)))

% clf; histc(rand(1000,1),0:0.1:1)
