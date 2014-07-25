% pModelSpecificCutoff(Coeff,Features,Method) evaluate sequence features in Features against linear criteria in Coeff

% Features has these entries:
% 1 AlignmentScore
% 2 CoreEditDistance
% 3 Percentile

function [Met,CutoffScore] = pModelSpecificCutoff(GroupData,Features,Params)

if nargin < 3,
	Params.DeficitCutoff    = 20;
	Params.CoreEditCutoff   = 5;
	Params.PercentileCutoff = 0.2;
	Params.CutoffType       = 2;
end

Met = 1;

switch Params.CutoffType,
case 1                                        % no cutoff

	Met = ones(size(Features(:,1)));

	if isfield(GroupData,'DeficitEditCutoff'),
		MixedScore = GroupData.CoreEditCoeff * Features(:,2) + GroupData.DeficitCoeff * (max(GroupData.OwnScore) - Features(:,1));
		CutoffScore =  100*(MixedScore - GroupData.DeficitEditCutoff)/(-GroupData.DeficitEditCutoff);
		CutoffScore = CutoffScore .* ((CutoffScore < 0) + (CutoffScore > 0) .* (Features(:,1) >= GroupData.MinScore) .* (Features(:,2) <= GroupData.CoreEditCutoff));
	else	
		CutoffScore = zeros(size(Met));
	end

case 2                                        % generic cutoffs only

%	Met = ((max(GroupData.OwnScore)-Features(:,1) <= Params.DeficitCutoff) .* (Features(:,2) <= Params.CoreEditCutoff) .* (Features(:,3) >= Params.PercentileCutoff));
	Met = ((max(GroupData.OwnScore)-Features(:,1) <= Params.DeficitCutoff) .* (Features(:,2) <= Params.CoreEditCutoff));
	Met = Met + (Features(:,2) == 0);           % core edit distance 0 means it's a match
	Met = min(1,Met);

	if isfield(GroupData,'DeficitEditCutoff'),
		MixedScore = GroupData.CoreEditCoeff * Features(:,2) + GroupData.DeficitCoeff * (max(GroupData.OwnScore) - Features(:,1));
		CutoffScore =  100*(MixedScore - GroupData.DeficitEditCutoff)/(-GroupData.DeficitEditCutoff);
		CutoffScore = CutoffScore .* ((CutoffScore < 0) + (CutoffScore > 0) .* (Features(:,1) >= GroupData.MinScore) .* (Features(:,2) <= GroupData.CoreEditCutoff));
	else	
		CutoffScore = zeros(size(Met));
	end

case 3                                        % model-specific cutoff

	MixedScore = GroupData.CoreEditCoeff * Features(:,2) + GroupData.DeficitCoeff * (max(GroupData.OwnScore) - Features(:,1));
	CutoffScore =  100*(MixedScore - GroupData.DeficitEditCutoff)/(-GroupData.DeficitEditCutoff);
	CutoffScore = CutoffScore .* ((CutoffScore < 0) + (CutoffScore > 0) .* (Features(:,1) >= GroupData.MinScore) .* (Features(:,2) <= GroupData.CoreEditCutoff));
	CutoffScore = max(-999,CutoffScore);                      % no need to go further than this

	Met = (CutoffScore > 0) + (Features(:,2) == 0);           % core edit distance 0 means it's a match
	Met = min(1,Met);

end

