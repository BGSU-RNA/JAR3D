% pIndividualGroupDiagnostic goes through each sequence,
% looks up the correct group and its cutoff score,
% sees if the sequence meets the cutoff score of its own group,
% and counts the number of other groups with a better score

function [NumBetterScore] = pIndividualGroupCutoffScore(GroupData,OnlyStructured,OwnMotif,Score,FASTA)

if OnlyStructured == 1
  str = find((cat(1,GroupData.Structured) == 1));  % structured internal loops
else
  str = 1:length(GroupData);
end

NumBetterScore = zeros(1,length(OwnMotif));

NumDontMatchOwnGroup = 0;

for n = 1:length(OwnMotif)                  % loop through sequences
  g = OwnMotif(n);                          % group that this sequence belongs to

  if Score(n,g,1) <= 0
    fprintf('Sequence %4d in group %s has cutoff score %10.4f; %s\n', n, GroupData(g).MotifID, Score(n,g,1), FASTA(n).Sequence);
    NumDontMatchOwnGroup = NumDontMatchOwnGroup + 1;
  end

  if OnlyStructured == 0 || GroupData(g).Structured > 0
    m = max(Score(n,:,:),[],3);             % maximum over all rotations of score against each group
    GroupsWithBetterCutoff = find((m(str) > Score(n,g,1)) .* (m(str) > 0) );
    NumBetterScore(n) = length(GroupsWithBetterCutoff);
    if NumBetterScore(n) > 0
      for k = GroupsWithBetterCutoff        % loop through pointers to groups with better score
        fprintf('Sequence %4d in group %s with cutoff score %10.4f does better against group %s cutoff score %10.4f; %s\n', n, GroupData(g).MotifID, Score(n,g,1), GroupData(str(k)).MotifID, m(str(k)), FASTA(n).Sequence);
      end
    end
  else
    NumBetterScore(n) = -1;                 % not structured
  end
end

