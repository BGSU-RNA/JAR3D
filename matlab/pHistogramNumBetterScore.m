
function [void] = pHistogramNumBetterScore(NumBetterScore,T,fs)

maxrank = 10;
n = hist(min(NumBetterScore(NumBetterScore >= 0),maxrank),0:maxrank);
if length(n) < 10,
  n(10) = 0;
end
bar(0:maxrank,n/sum(n));
axis([-0.5 maxrank+0.5 0 max(n/sum(n))*1.1]);
text(4,max(n/sum(n))*0.8,[sprintf('%0.2f',100*n(1)/sum(n)) '% had the correct model']);
text(4,max(n/sum(n))*0.7,[sprintf('%0.2f',100*n(2)/sum(n)) '% had 1 other model score better']);
text(4,max(n/sum(n))*0.6,[sprintf('%0.2f',100*n(3)/sum(n)) '% had 2 other models score better']);
text(4,max(n/sum(n))*0.5,[sprintf('%0.2f',100*n(4)/sum(n)) '% had 3 other models score better']);
text(4,max(n/sum(n))*0.4,[sprintf('%0.2f',100*n(5)/sum(n)) '% had 4 other models score better']);
set(gca,'XTick',0:10)
set(gca,'XTickLabel',{'0','1','2','3','4','5','6','7','8','9','>=10'})
title(T,'fontsize',fs);
xlabel(['Number of models that score higher than the correct model'],'fontsize',fs);
ylabel('Frequency','fontsize',fs);
