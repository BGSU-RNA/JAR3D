% pIndividualGroupDiagnostic goes through each sequence,
% looks up the correct group and its score,
% and counts the number of other groups with a better score

function [NumBetterScore] = pIndividualGroupDiagnostic(GroupData,OnlyStructured,OwnMotif,Score)

if OnlyStructured == 1
  str = find((cat(1,GroupData.Structured) == 1));  % structured internal loops
else
  str = 1:length(GroupData);
end

for n = 1:length(OwnMotif)                  % loop through sequences
  g = OwnMotif(n);                          % group that this sequence belongs to
  if OnlyStructured == 0 || GroupData(g).Structured > 0
    m = max(Score(n,:,:),[],3);             % maximum over all rotations of score against each group
    NumBetterScore(n) = length(find(m(str) > Score(n,g,1)));
    if NumBetterScore(n) > 0

      for k = find(m(str) > Score(n,g,1))    % loop through groups with better score
        fprintf('Sequence %d in group %s with score %10.6f does better against group %s score %10.6f\n', n, GroupData(g).MotifID, Score(n,g,1), GroupData(str(k)).MotifID, m(str(k)));
      end
    end
  else
    NumBetterScore(n) = -1;                 % not structured
  end
end

