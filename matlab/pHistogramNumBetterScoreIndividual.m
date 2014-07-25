
function [counts] = pHistogramNumBetterScoreIndividual(NumBetterScore,cc,T,fs,FASTA)

% if supplied, use multiplicity information for sequences

if nargin >= 5,
  NBS = [];
  newcc = [];
  for k = 1:length(NumBetterScore),
  	m = FASTA(k).Multiplicity;
    NBS = [NBS ones(1,m)*NumBetterScore(k)];
    newcc = [newcc ones(1,m)*cc(k)];
  end

  NumBetterScore = NBS;
  cc = newcc;  
end


% remove sequences with negative NumBetterScore (code to ignore)

j = find(NumBetterScore >= 0);
NumBetterScore = NumBetterScore(j);
cc = cc(j);


maxrank = 10;
n = hist(min(NumBetterScore,maxrank),0:maxrank);
if length(n) < 10,
  n(10) = 0;
end
%bar(0:maxrank,n/sum(n));
clf
counts = zStackedHistogram(min(NumBetterScore,maxrank),abs(cc),1);

text(4,max(n/sum(n))*0.8,[sprintf('%0.2f',100*n(1)/sum(n)) '% had the correct model']);
text(4,max(n/sum(n))*0.7,[sprintf('%0.2f',100*n(2)/sum(n)) '% had 1 other model score better']);
text(4,max(n/sum(n))*0.6,[sprintf('%0.2f',100*n(3)/sum(n)) '% had 2 other models score better']);
text(4,max(n/sum(n))*0.5,[sprintf('%0.2f',100*n(4)/sum(n)) '% had 3 other models score better']);
text(4,max(n/sum(n))*0.4,[sprintf('%0.2f',100*n(5)/sum(n)) '% had 4 other models score better']);
text(4,max(n/sum(n))*0.3,[sprintf('%0.2f',100*n(11)/sum(n)) '% had 10 or more other models score better']);

% Note:  bars are colored by absolute difference in sequence length
n = counts(:,end);
axis([-0.5 maxrank+0.5 0 max(n/sum(n))*1.05]);
set(gca,'XTick',0:10)
set(gca,'XTickLabel',{'0','1','2','3','4','5','6','7','8','9','>=10'})
title(T,'fontsize',fs);
xlabel(['Number of models that score higher than correct model'],'fontsize',fs);
