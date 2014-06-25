% pMultipleGroupDiagnostic randomly selects 1, 2, 3, etc. sequences from
% each sequence group and sees how many motif groups score better than
% the sequence group's own group

% MLPS(a,m,r) is the score of sequence a against model m using rotation r

NumSeqsEach = 1;                              % number of sequences to average
NumSeqsEach = max(GroupSize);
MaxNumReps = max(GroupSize);
MaxNumReps = 100;
MaxNumReps = 1;

NSEvalues = [1:10 15 20 25 30 35 40 45 50 max(GroupSize)];

fff = 5*ones(size(NSEvalues));
fff(end) = 1;

for bbb = 1:length(NSEvalues),

NSE = NSEvalues(bbb);

clear NumBetterScore
n = 1;

for g = 1:length(GroupData),                          % loop through groups
  j = find(OwnMotif == g);                    % all sequences in group n

  % reproduce the group-group diagnostic, weighing each group the same
  NumSeqsEach = max(GroupSize);               % use all sequences each time
  NumReps     = 1;

  % reproduce the individual-group diagnostic weighing each group the same
  NumSeqsEach = 1;
  NumReps     = max(GroupSize);               % 

  % nearly reproduce the group-group diagnostic, weighing by sequence
  NumSeqsEach = max(GroupSize);               % use all sequences each time
  NumReps     = length(j);

  % nearly reproduce the individual-group diagnostic, weighing by sequence
  NumSeqsEach = NSE;
  NumReps     = 25*length(j);

  % weigh each group the same
  NumSeqsEach = NSE;                          % different numbers of seqs
  NumReps     = fff(bbb)*ceil(length(j)/NSE); % sliding scale

  for r = 1:NumReps,
    mlps = zeros(size(MLPS(1,:,:)));          % place to store scores
    p = randperm(length(j));                  % randomly order sequences
    p = p(1:min(length(j),NumSeqsEach));      % use first few sequences
    for k = 1:length(p),                      % loop through chosen sequences
      mlps = mlps + MLPS(j(p(k)),:,:);        % sum scores of all sequences
    end
    m = max(mlps,[],3);                       % look for biggest sum over rotations
    NumBetterScore(n) = length(find(m > mlps(1,g,1)));
    n = n + 1;
  end
end

NBS = NumBetterScore;

% fprintf('%d comparisons were made\n', length(NBS));

figure(17)
clf
maxrank = 10;
n = hist(min(NBS(NBS >= 0),maxrank),0:maxrank);
if length(n) < 10,
  n(10) = 0;
end
bar(0:maxrank,n/sum(n));
axis([-0.5 maxrank+0.5 0 max(n/sum(n))*1.1]);
set(gca,'XTick',0:10)
set(gca,'XTickLabel',{'0','1','2','3','4','5','6','7','8','9','>=10'})
if NSE == 1,
  text(4,max(n/sum(n))*0.9,['Using ' num2str(NSE) ' sequence at a time,']);
else
  text(4,max(n/sum(n))*0.9,['Using ' num2str(NSE) ' sequences at a time,']);
end
text(4,max(n/sum(n))*0.8,[sprintf('%0.1f',100*n(1)/sum(n)) '% had the correct model score highest']);
text(4,max(n/sum(n))*0.7,[sprintf('%0.1f',100*n(2)/sum(n)) '% had 1 other model score better']);
text(4,max(n/sum(n))*0.6,[sprintf('%0.1f',100*n(3)/sum(n)) '% had 2 other models score better']);
text(4,max(n/sum(n))*0.5,[sprintf('%0.1f',100*n(4)/sum(n)) '% had 3 other models score better']);
text(4,max(n/sum(n))*0.4,[sprintf('%0.1f',100*n(5)/sum(n)) '% had 4 other models score better']);
text(4,max(n/sum(n))*0.3,[sprintf('%0.1f',100*n(11)/sum(n)) '% had 10 or more other models score better']);
title(['Multiple-group diagnostic, ' num2str(length(SeqNames)) ' sequence groups against ' num2str(length(GroupData)) ' models'],'fontsize',fs);
xlabel(['Number of models that score higher than the correct model'],'fontsize',fs);
ylabel('Frequency','fontsize',fs);
if NSE == max(GroupSize),
  print(gcf,'-dpng',[DiagnosticPath filesep 'Multiple_Group_Max.png']);
else
  print(gcf,'-dpng',[DiagnosticPath filesep 'Multiple_Group_' num2str(NumSeqsEach) '.png']);
end

end

% old code to determine the number of repetitions

if 0 > 1,
  if length(j) <= NumSeqsEach,
    NumReps = 1;
  else
    NumReps = MaxNumReps;
  end

  if NumSeqsEach == 1,
    NumReps = length(j);                      % reproduces sequence-weighted individual-group diagnostic
  end
end
