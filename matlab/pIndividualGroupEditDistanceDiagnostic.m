
function [NumBetterScore] = pIndividualGroupEditDistanceDiagnostic(Structured,OnlyStructured,OwnMotif,EditDistanceIndividual)

if OnlyStructured == 1,
  str = find((Structured == 1));  % structured internal loops
else
  str = 1:length(Structured);
end

for n = 1:length(OwnMotif),                 % loop through sequences
  g = OwnMotif(n);                          % correct motif group
  if OnlyStructured == 0 || Structured(g) > 0,
    EDI = EditDistanceIndividual(n,:,:);    % extract distances for this seq
    m = min(EDI,[],3);                      % best over rotations
    if any(EditDistanceIndividual(n,g,:) == Inf),
      NumBetterScore(n) = -1;               % exact sequence match in one rotation
                                            % no different core sequence!
    else
      NumBetterScore(n) = length(find(m(str) < EditDistanceIndividual(n,g,1)));
    end
  else
    NumBetterScore(n) = -1;                 % not structured
  end
end

