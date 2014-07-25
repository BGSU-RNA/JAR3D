
function [NBS] = pGroupGroupDiagnostic(NumSeqGroups,GroupData,OnlyStructured,FASTA,MLPS,OwnMotif,SeqGroup)

NumBetterScore = zeros(1,NumSeqGroups);

if OnlyStructured == 1,
  str = find(cat(1,GroupData.Structured) == 1);  % structured loops
else
  str = 1:length(GroupData);
end

for n = 1:NumSeqGroups,                       % loop through groups
  j = find(SeqGroup == n);                    % all sequences in group n
  if ~isempty(j),
    if OnlyStructured == 0 || GroupData(OwnMotif(j(1))).Structured == 1,
      mlps = zeros(size(MLPS(1,:,:)));
      for m = 1:length(j),
        mult = FASTA(j(m)).Multiplicity;
        mlps = mlps + mult * MLPS(j(m),:,:);  % sum scores of all sequences
      end
      m = max(mlps,[],3);                     % look for biggest sum
      g = OwnMotif(j(1));
      NumBetterScore(n) = length(find(m(str) > mlps(1,g,1)));
    else
      NumBetterScore(n) = -1;                 % group is not structured
    end
  else
    NumBetterScore(n) = -1;                   % group is not represented
  end
end

NBS = NumBetterScore(NumBetterScore >= 0);
