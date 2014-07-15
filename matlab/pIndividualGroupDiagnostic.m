% pIndividualGroupDiagnostic goes through each sequence, looks up the correct group and its score, and counts the number of other groups with a better score

function [NumBetterScore] = pIndividualGroupDiagnostic(GroupData,OnlyStructured,OwnMotif,Score)

if OnlyStructured == 1,
  str = find((cat(1,GroupData.Structured) == 1));  % structured internal loops
else
  str = 1:length(GroupData);
end

for n = 1:length(OwnMotif),                     % loop through sequences
  g = OwnMotif(n);
  if OnlyStructured == 0 || GroupData(g).Structured > 0,
    m = max(Score(n,:,:),[],3);                 % maximum over rotations, if any
    NumBetterScore(n) = length(find(m(str) > Score(n,g,1)));
  else
    NumBetterScore(n) = -1;                 % not structured
  end
end

