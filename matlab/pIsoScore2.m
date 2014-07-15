% pIsoScore2 returns a 4x4 matrix of scores depending on the specified
% interaction category Class and the specified Pair having this interaction.

function [S,IDI] = pIsoScore2(Class,Code1,Code2,ExemplarIDI,Normalize)

if nargin < 4,
  load PairExemplars
end

if strcmp(class(Class),'char'),
  Class = xGetEdgeNums(Class);
end

if strcmp(class(Code1),'char'),
  Code1 = pLtoN(Code1);
end

if strcmp(class(Code2),'char'),
  Code2 = pLtoN(Code2);
end

% ---------------------------------------- look up IDI values

% [Class Code1 Code2]

IDI = zExemplarIDI(Class,Code1,Code2,ExemplarIDI);

% ---------------------------------------- turn IDI into scores

for a = 1:4,
  for b = 1:4,
    if isnan(IDI(a,b)),
      IDI(a,b) = 100;                    % this pair is not possible
    end
  end
end

%S = 1 * (IDI < 2.0) + 1 * exp(-0.2*(IDI-2.0)) .* (IDI >= 2.0);

% -------------- Map IDI to scorue using smooth quadratic-like function

S = 1 ./ (1 + (0.5*IDI).^2);                   % score decreases with IDI
S = max(S,0.05);                         % set a minimum score for all pairs

% -------------- Map IDI to score using piecewise linear function

 x = [0 1.8 2.20 3.10 3.50 7.00 9.00 Inf];
 y = [1 0.8 0.13 0.07 0.05 0.03 0.01 0.01];
%y = [1 0.8 0.30 0.15 0.05 0.03 0.01 0.01];

if 0 > 1,
  x = [0 1.8 2.20 3.10 3.50 7.00 9.00 20 Inf];
  y = [1 0.8 0.13 0.07 0.05 0.03 0.01 0.01 0.01];

  figure(1)
  clf
  plot(x,y,'k','linewidth',2);  
  fs = 18;
  xlabel('IsoDiscrepancy Index (IDI)','fontsize',fs)
  ylabel('Probability score','fontsize',fs);
  set(gca,'fontsize',fs)
  axis([0 10 0 1])
  saveas(gcf,'C:\Users\zirbel\Dropbox\2014 JAR3D article\IDI to Score mapping.png')
end

for a = 1:4,
  for b = 1:4,
    w = find(IDI(a,b) < x);
    w = w(1);
    S(a,b) = y(w-1) + (IDI(a,b)-x(w-1))*(y(w)-y(w-1))/(x(w)-x(w-1));
  end
end

% ---------------------------------------- normalize scores into probabilities

if Normalize == 1
    S = S / sum(sum(S));                  % normalize
end

% -------------- Code to test mapping IDI to score using piecewise linear

if 0 > 1,
  clear R
  clear S
  clf
  for v = 0:1000,
    IDI = 10*v/1000;
    a = find(IDI < x);
    a = a(1);
    S(v+1) = y(a-1) + (IDI-x(a-1))*(y(a)-y(a-1))/(x(a)-x(a-1));
    R(v+1) = IDI;
  end
  plot(R,S);  
  hold on
  IDI = (0:1000)/100;
  S = 1 ./ (1 + (0.5*IDI).^2);                   % score decreases with IDI
  S = max(S,0.05);                         % set a minimum score for all pairs

  plot(IDI,S,'r');  
  hold on
  xlabel('IDI value, 2-3.3 is near isosteric');
  ylabel('Probability score')
end
