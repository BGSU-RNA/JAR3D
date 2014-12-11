% pSetModelSpecificCutoffs sets model-specific cutoffs using data from alignments and randomly-generated data

% Each sequence is represented by a structured variable whose fields are:
% .Deficit
% .CoreEdit
% .FullEdit
% .Conserved
% .NumInstances
% .MLPS
% .Length
% .Source    = 1 for alignment, = 2 for random but which maximizes alignment score, 2.2 for worse random sequences, etc.
% .Sequence
% .LoopID
% .MotifID
% .NumFixed
% .Basepairs

function [void] = pSetModelSpecificCutoffs(OutputBase,Release,UseAlignmentSequences)

if nargin < 3,
	UseAlignmentSequences = 2;               % use alignment sequences to set cutoffs and show them too
	UseAlignmentSequences = 1;               % show alignment sequences but don't use them to set cutoffs
	UseAlignmentSequences = 0;               % don't use or show alignment sequences
end

Params.DeficitCutoff    = 20;              % don't ever change this number because it has become hard-coded into the java code
Params.CoreEditCutoff   = 5;

Grayscale = 0;                              % plot in grayscale for the paper or in color for talks
tfs = 13;

if UseAlignmentSequences == 2,
	NumAlignmentSequencesNeeded = 20;           % number for cutoffs based on alignment and random sequences
else
	NumAlignmentSequencesNeeded = Inf;          % for cutoffs based only on random sequences
end

Show2dPlots = 0;
Show3dPlots = 0;
Save3DPlots = 0;
SaveDeficitCoreEditPlot = 1;
SaveGroupData = 1;
ListSequences = 0;
PlotCutoffPlane = 0;
CoreEditCutoff = 5;
MaxNumFP = [200 200 400 500 600 800 1000];  % maximum number of false positives according to number of conserved basepairs
MaxNumFP = [MaxNumFP Inf * ones(1,200)];
Params.CutoffType = 2;

if Save3DPlots > 0,
	Show3dPlots = 1;
end

Release = strrep(Release,'/',filesep);
Release = strrep(Release,'\',filesep);
OutputPath = [OutputBase filesep Release];
MotifRelease = Release(4:end);

loopType = Release(1:2);

ModelPath = [OutputPath filesep 'lib'];

GroupData = pGetModelData(OutputPath,loopType);

MSCOutputPath = [OutputPath filesep 'ModelSpecificCutoffs'];

if ~(exist(MSCOutputPath) == 7),        % if directory doesn't yet exist
  mkdir(MSCOutputPath);
end

DiaryFile = [MSCOutputPath filesep 'log_' date '.txt'];

diary(DiaryFile);

if ~(exist(MSCOutputPath) == 7),        % if directory doesn't yet exist
  mkdir(MSCOutputPath);
end

if UseAlignmentSequences > 0 && ~exist('AlignmentData'),
	load([OutputPath filesep loopType '_Alignment_Sequence_Data.mat']);
	AlignmentData = SequenceData;
	clear SequenceData

	if isfield(AlignmentData,'Percentile'),
		AlignmentData = rmfield(AlignmentData,'Percentile');
	end

	AlignmentDataMotifIDs = cell(1,length(AlignmentData));
	for i = 1:length(AlignmentData),
		AlignmentDataMotifIDs{i} = AlignmentData(i).MotifID;
	end

	fprintf('Loaded data from %d sequences from alignments\n',length(AlignmentData));
else
	AlignmentData = [];
	AlignmentDataMotifIDs = {};
end

[y,motiforder] = sort(cat(1,GroupData.NumNT),1,'ascend');

motiforder = 1:length(GroupData);

for iii = 1:length(GroupData),

	motifnum = motiforder(iii);
	pauseafter = 0;

	if ListSequences > 0,
		clc
	end

	GroupData(motifnum)
	
	CurrentMotif = GroupData(motifnum).MotifID;

	% -------------------------- filter out poor alignment data

	k = find(ismember(AlignmentDataMotifIDs,CurrentMotif));
  SequenceData = AlignmentData(k);

  if length(SequenceData) > 0,
		ce = cat(1,SequenceData.CoreEdit);
		k = find(ce <= CoreEditCutoff);
		SequenceData = SequenceData(k);

		def = cat(1,SequenceData.Deficit);
		k = find(def <= Params.DeficitCutoff);
		SequenceData = SequenceData(k);
	  fprintf('Removed based on deficit, %d sequences remaining\n',length(SequenceData));
		fprintf('%10d sequences from alignments with deficits between %d %d\n', length(SequenceData), 0, Params.DeficitCutoff);

		[y,k] = sort(rand(1,length(SequenceData)));         % randomize order
		SequenceData = SequenceData(k);
	end

	% ---------------------------- record data from motif group

	clear NSD
	for co = 1:GroupData(motifnum).NumInstances,
		NSD(co).Deficit      = max(GroupData(motifnum).OwnScore) - GroupData(motifnum).OwnScore(co);
		NSD(co).CoreEdit     = 0;
		NSD(co).FullEdit     = 0;
		NSD(co).Conserved    = GroupData(motifnum).NumNT;
		NSD(co).NumInstances = GroupData(motifnum).NumInstances;
		NSD(co).MLPS         = GroupData(motifnum).OwnScore(co);
		NSD(co).Length       = length(GroupData(motifnum).OwnSequence{co}) - length(strfind(GroupData(motifnum).OwnSequence{co},'*'));
		NSD(co).Source       = 0;                                    % from 3D structures
		NSD(co).Basepairs    = GroupData(motifnum).NumBasepairs;
		NSD(co).NumFixed     = GroupData(motifnum).NumFixed;
		NSD(co).MotifID      = GroupData(motifnum).MotifID;
		NSD(co).Sequence     = GroupData(motifnum).OwnSequence{1};
	end

	fprintf('%d sequences from 3D structures\n',length(NSD));

	% ---------------------------- record data from random sequences

	AllSD = [NSD SequenceData];      % put sequences from 3D first

	if isfield(GroupData(motifnum),'DeficitEditData'),
		RSNum = length(GroupData(motifnum).DeficitEditData(:,1));
	else
		fprintf('pSetModelSpecificCutoffs:  No data about how random sequences match the model\n');

		for mm = 1:length(GroupData),
			if ~isfield(GroupData(mm),'DeficitEditData'),
				fprintf('Group %3d is missing DeficitEditData\n',mm);
			end
		end

		pause
		RSNum = 0;
	end

	fprintf('Using %d random sequences, %d from an alignment, and %d from 3D structures\n', RSNum, length(SequenceData), length(NSD));

	if RSNum + length(AllSD) > 0,

		clear SequenceData

		Keep = ones(1,length(AllSD));
		for i = 1:length(AllSD),
			if isempty(AllSD(i).CoreEdit),
				Keep(1,i) = 0;
			end
		end
		AllSD = AllSD(find(Keep));

		Source = cat(1,AllSD.Source);

		clear SD
		SD(:,1) = cat(1,AllSD.Deficit);
		SD(:,2) = cat(1,AllSD.CoreEdit);
		SD(:,3) = cat(1,AllSD.FullEdit);
		SD(:,4) = cat(1,AllSD.Length);
		SD(:,5) = cat(1,AllSD.Conserved);
		SD(:,6) = SD(:,3) - SD(:,2);
		SD(:,7) = cat(1,AllSD.MLPS);
		SD(:,10) = cat(1,AllSD.NumFixed);
		SD(:,11) = cat(1,AllSD.Basepairs);

		if RSNum > 0,
			clear RSD
			RSD(:,1) = cat(1,GroupData(motifnum).DeficitEditData(:,1));
			RSD(:,2) = cat(1,GroupData(motifnum).DeficitEditData(:,2));
			RSD(:,3) = 99 * ones(RSNum,1);                                % just pretend
			RSD(:,4) = zeros(RSNum,1);                                    % just pretend
			RSD(:,5) = GroupData(motifnum).NumNT * ones(RSNum,1);
			RSD(:,6) = RSD(:,3) - RSD(:,2);
			RSD(:,7) = max(GroupData(motifnum).OwnScore) - RSD(:,1);
			RSD(:,10) = GroupData(motifnum).NumFixed * ones(RSNum,1);
			RSD(:,11) = GroupData(motifnum).NumBasepairs * ones(RSNum,1);

			SD = [SD; RSD];              % append randomly-generated sequences
		end

		Source = [Source; 2*ones(RSNum,1)];

		LL = length(Source);

		clear SDR
		SDR(:,1) = SD(:,1);
		SDR(:,2) = SD(:,2) + 0.25*max(0,(min(Source,2)-1)) + 0.25*rand(LL,1);
		i = find(Source == 0);
		SDR(i,2) = SD(i,2);                  % no random shift for edit distances for sequences from 3D structures
		SDR(:,3) = SD(:,3) + rand(LL,1)/4;
		SDR(:,4) = SD(:,4) + rand(LL,1)/4;
		SDR(:,5) = SD(:,5) + rand(LL,1)/4;
		SDR(:,6) = SD(:,6) + rand(LL,1)/4;
		SDR(:,7) = SD(:,7);
		SDR(:,10) = SD(:,10) + rand(LL,1)/4;
		SDR(:,11) = SD(:,11) + rand(LL,1)/4;

		k = find((SD(:,2) == 0) .* (Source==2));        % random sequences with core edit distance 0
		Source(k) = 3;

		Names{1} = 'Alignment Score Deficit';
		Names{2} = 'Interior Edit Distance';
		Names{3} = 'FullEdit';
		Names{4} = 'SequenceLength';
		Names{5} = 'NumConserved';
		Names{6} = 'FlankEdit';
		Names{7} = 'AlignmentScore';
	  Names{10} = 'NumFixed';
	  Names{11} = 'Basepairs';

		f = 1;
	  if Show2dPlots > 0,
			j = 1:min(length(Source),10000);

		  n = [1 2 3 4 6 7];                  % which variables to plot and in what order

			for a = 1:length(n),
				for b = (a+1):length(n),
					figure(f)
					clf
					scatter(SDR(j,n(a)),SDR(j,n(b)),4,min(Source(j),3),'filled');
					caxis([1 4])
					xlabel(Names{n(a)});
					ylabel(Names{n(b)});
					f = f + 1;
				end
			end
		end

		if Show3dPlots > 0,

			Triple(1,:) = [1 2 8];
			Triple(2,:) = [1 2 9];
			Triple(3,:) = [1 2 4];
			Triple(4,:) = [4 5 1];
			Triple(5,:) = [4 5 7];
			Triple(6,:) = [3 5 11];
			Triple(7,:) = [4 5 10];
			Triple(8,:) = [2 4 7];

	  	clear Triple

			Triple(1,:) = [2 4 1];

			j = 1:min(2000,length(Source));  % plot up to 2000 points
			j = [j find(Source == 0)'];

			for t = 1:length(Triple(:,1)),
				figure(f)
				[ax,az] = view;
				clf
				clf
				scatter3(SDR(j,Triple(t,1)),SDR(j,Triple(t,2)),SDR(j,Triple(t,3)),4+40*(Source(j) == 0),min(Source(j),2),'filled')
				xlabel(Names{Triple(t,1)});
				ylabel(Names{Triple(t,2)});
				zlabel(Names{Triple(t,3)});
				rotate3d on
				caxis([-0.5 3]);
				view([ax az]);                % keep the current orientation
				drawnow
				f = f + 1;
				axis([0 5.2 0 20 0 20]);                           % show every case on the same scale
				view([0 0]);
			end
		end

		DeficitCoeff = GroupData(motifnum).DeficitCoeff;
		GroupData(motifnum).CoreEditCoeff = 3;               % insist on this being 3
		CoreEditCoeff = GroupData(motifnum).CoreEditCoeff;
		Coeff = [DeficitCoeff CoreEditCoeff];

		BinarySource = 1*(Source <= 1) + 2*(Source > 1);

		i = [2 4 1];  % core edit, sequence length, deficit

  	if ListSequences > 0,
  		OrderMatrix = [BinarySource SD(:,7)];
  		[y,ii] = sortrows(OrderMatrix,[1 2]);
  		for n = 1:length(ii),
  			j = ii(n);
  			fprintf('Source %d %50s CoreEdit %d Score %8.4f Deficit %8.4f \n', Source(j), AllSD(j).Sequence, SD(j,2), SD(j,7), SD(j,1));
  		end
  	end

		ra = find(Source == 1);
    AlignmentMixedScores = SD(ra,[1 2]) * Coeff';          % scores of these sequences against current coefficients

		ra = find(Source == 2);
    RandomMixedScores = SD(ra,[1 2]) * Coeff';             % scores of these sequences against current coefficients

		ra = find((Source == 2) .* (SD(:,2) > 1));             % random sequences at core edit distance 2 or greater
    RandomMixedScoresB = SD(ra,[1 2]) * Coeff';            % scores of these sequences against current coefficients

		GroupData(motifnum).RandomQuantile80 = quantile(RandomMixedScores,0.2);
		GroupData(motifnum).RandomQuantile96 = quantile(RandomMixedScores,0.04);
		GroupData(motifnum).RandomQuantile98 = quantile(RandomMixedScores,0.02);

		switch loopType,
		case 'IL'
			DefaultMixedScoreCutoff = min(25,10+1.8*(GroupData(motifnum).NumNT-5));
		case 'HL'
			DefaultMixedScoreCutoff = min(25,10+3*(GroupData(motifnum).NumNT-5));
		end
		MixedScoreCutoff = DefaultMixedScoreCutoff;
 		GroupData(motifnum).CutoffMethod = 1;

    figure(2)
  	clf

  	k = find((Source <= 1) + (Source == 2).*(SD(:,2) > 0));  % only use alignment sequences and "bad" random sequences
  	MixedScores = SD(k,[1 2]) * Coeff';

		if sum(Source == 2) > 0 && sum(Source <= 1) > NumAlignmentSequencesNeeded,  % enough of both random and alignment sequences

			% find the cutoff which maximizes the TN - (1-TP) rate, if large enough
	  	MixedScoreValues = 0:0.1:40;
	  	performance = zeros(length(MixedScoreValues),3);

	  	for c = 1:length(MixedScoreValues),
	  	  mynewclasses = 2 - (MixedScores <= MixedScoreValues(c));
	  	  TP = sum((mynewclasses == 1) .* (BinarySource(k,:) == 1)) / sum(BinarySource(k,:) == 1);
	  	  TN = sum((mynewclasses == 2) .* (BinarySource(k,:) == 2)) / sum(BinarySource(k,:) == 2);
				performance(c,:) = [TP TN TN-(1-TP)];
	  	end

	  	[mm,ii] = max(performance(:,3));

			if mm > 0.1 && mm < 0.4
				fprintf('Borderline case for imposing a model-specific cutoff, maximum %8.4f\n',mm);
			end

			if min(performance(:,3)) < -0.5,
				fprintf('Many random sequences perform better than sequences from alignments.\n');
			end

	  	figure(2)
	  	clf
	  	plot(MixedScoreValues,performance(:,3));
	  	hold on
	  	plot(MixedScoreValues,mm-0.02,'r');
	  	xlabel('MixedScore before adjustments')
	  	drawnow

	  	if mm > 0.2,                                              % random matches occur where true positives do not
	  		ii = sort(find(performance(:,3) > mm-0.02));            % find places with nearly as good a performance
	  		MixedScoreCutoff = MixedScoreValues(ii(1));                          % use the most restrictive one
	  		fprintf('MixedScore cutoff set to maximize TN - (1 - TP)\n');
	  		fprintf('TP rate %8.4f, TN rate %8.4f\n', 100*performance(ii(1),1), 100*performance(ii(1),2));
	  		GroupData(motifnum).CutoffMethod = 2;
	  	else
		    MixedScoreCutoff = (quantile(RandomMixedScores,0.04)+quantile(RandomMixedScoresB,0.02))/2;  % find the cutoff which gets the true negative rate roughly equal to 98%
	  		fprintf('Setting cutoff by percentiles of random sequences, %8.4f\n',MixedScoreCutoff);
	  		GroupData(motifnum).CutoffMethod = 3;
	  	end

	  	LocalMaxNumFP = MaxNumFP(GroupData(motifnum).NumBasepairs);

	  	NumFalsePositive = sum((2 - (MixedScores <= MixedScoreCutoff) == 1) .* (BinarySource(k,:) == 2));
 		  mynewclasses = 2 - (MixedScores <= MixedScoreCutoff);
		  TP = sum((mynewclasses == 1) .* (BinarySource(k,:) == 1)) / sum(BinarySource(k,:) == 1);

	  	if NumFalsePositive > LocalMaxNumFP && TP > 0.5,  % large number of false positives and an acceptable true positive rate
	  		while NumFalsePositive >= LocalMaxNumFP && MixedScoreCutoff > 0,
	  			fprintf('Constant %8.4f False positives %4d\n',MixedScoreCutoff,sum((2 - (MixedScores <= MixedScoreCutoff) == 1) .* (BinarySource(k,:) == 2)));
		  		MixedScoreCutoff = MixedScoreCutoff - 0.1;
			  	NumFalsePositive = sum((2 - (MixedScores <= MixedScoreCutoff) == 1) .* (BinarySource(k,:) == 2));
	  		end
    		fprintf('Tightened the cutoff to reduce the sheer number of false positives\n');
	  		GroupData(motifnum).CutoffMethod = GroupData(motifnum).CutoffMethod + 2;
	  	end

	  elseif sum(Source == 2) > 20,                           % 3D and random sequences, few alignment sequences
	    MixedScoreCutoff = (quantile(RandomMixedScores,0.04)+quantile(RandomMixedScoresB,0.02))/2;  % find the cutoff which gets the true negative rate roughly equal to 98%
	    MixedScoreCutoff = quantile(RandomMixedScores,0.05);  % find the cutoff which gets the true negative rate roughly equal to 96%
	    MixedScoreCutoff = quantile(RandomMixedScores,0.04);  % find the cutoff which gets the true negative rate roughly equal to 96%
  		GroupData(motifnum).CutoffMethod = 6;

	    if MixedScoreCutoff > 20,
				fprintf('Decreased cutoff from %8.4f because the cutoff seemed overly generous\n',MixedScoreCutoff);
		    MixedScoreCutoff = max(20,min(RandomMixedScores)-0.1);
		    MixedScoreCutoff = max(20,quantile(RandomMixedScores,0.02));  % find the cutoff which gets the true negative rate roughly equal to 96%
	  		GroupData(motifnum).CutoffMethod = 8;
			end
 
	  elseif sum(Source == 2) == 0,                         % no random sequences meet the generic cutoffs
	  	MixedScoreCutoff = 35;                                        % impose essentially no model-specific cutoff
	  	fprintf('No random sequence matches, so essentially no model-specific cutoff imposed\n');
  		GroupData(motifnum).CutoffMethod = 9;

	  elseif sum(Source == 1) > NumAlignmentSequencesNeeded,                       % 3D and alignment sequences, few random sequences
	    MixedScoreCutoff = quantile(AlignmentMixedScores,0.90);                    % 90th percentile of alignment sequences
	  	GroupData(motifnum).CutoffMethod = 10;
		end

		if MixedScoreCutoff < 9.5,
			MixedScoreCutoff = 5 + Coeff(2);
			MixedScoreCutoff = 9.5;
			fprintf('Increased cutoff to %8.4f so that at least a few sequences with core edit distance 1 can meet the cutoff\n',MixedScoreCutoff);
	  	GroupData(motifnum).CutoffMethod = 11;
		end

		if MixedScoreCutoff > 25,
			MixedScoreCutoff = 25;
			fprintf('Decreased cutoff to %8.4f so that it is possible to reject matches\n',MixedScoreCutoff);
	  	GroupData(motifnum).CutoffMethod = 12;
		end

	  mynewclasses = 2 - (MixedScores <= MixedScoreCutoff);
	  TP = sum((mynewclasses == 1) .* (BinarySource(k,:) == 1)) / sum(BinarySource(k,:) == 1);
	  TN = sum((mynewclasses == 2) .* (BinarySource(k,:) == 2)) / sum(BinarySource(k,:) == 2);

		GroupData(motifnum).TruePositiveRate           = TP;
		GroupData(motifnum).TrueNegativeRate           = TN;
		GroupData(motifnum).NumberOf3DSequences        = length(GroupData(motifnum).OwnScore);
		GroupData(motifnum).NumberOfAlignmentSequences = sum(Source == 1);
		GroupData(motifnum).NumberOfRandomSequences    = sum(Source == 2);
		GroupData(motifnum).NumberOfFalsePositives     = sum((mynewclasses == 1) .* (Source(k,:) == 2));  % number of false positives

		GroupData(motifnum).DeficitCutoff    = Params.DeficitCutoff;
		GroupData(motifnum).CoreEditCutoff   = Params.CoreEditCutoff;

		GroupData(motifnum).DeficitEditCutoff = MixedScoreCutoff;

		GroupData(motifnum).MinScore        = max(GroupData(motifnum).OwnScore) - Params.DeficitCutoff;
		GroupData(motifnum).MaxScore        = max(GroupData(motifnum).OwnScore);
		GroupData(motifnum).ScoreEditCutoff = GroupData(motifnum).DeficitEditCutoff - Coeff(1) * max(GroupData(motifnum).OwnScore);

		m = motifnum;
		GD = GroupData(motifnum);
		fprintf('Group %3d, %-11s has acceptance rules AlignmentScore >= %8.4f, CoreEdit <= %d, and %8.4f * CoreEdit - %8.4f * AlignmentScore <= %8.4f\n',motifnum, GD.MotifID, GD.MinScore, GD.CoreEditCutoff, GD.CoreEditCoeff, GD.DeficitCoeff, GD.ScoreEditCutoff);
		fprintf('TP %8.2f%%, TN %8.2f%%, min %8.2f%%, %3d 3D sequences, %5d alignment sequences, %4d random sequences, %4d random matches, %2d NTs, %s\n',100*GD.TruePositiveRate, 100*GD.TrueNegativeRate, 100*min(GD.TruePositiveRate,GD.TrueNegativeRate), GD.NumInstances, GD.NumberOfAlignmentSequences, GD.NumberOfRandomSequences, GD.NumberOfFalsePositives, GD.NumNT, GD.Signature{1});

    % ---------------------------- Deficit versus core edit distance - the graph that really matters
		% ---------------- main figure - represent points from 3D, from alignments, random sequences
    figure(1)
    clf

    i = find(Source >= 2);                        % randomly-generated sequences
    if length(i) > 2000,
    	i = i(1:2000);
    end
    ii = find(Source == 1);                        % sequences from alignments
    if length(i) > 2000,
    	ii = ii(1:2000);
    end
    iii = find(Source == 0);                        % use all sequences from 3D

		ced = 0.01:0.01:100;                                % range of core edit distances
  	ycutoff = (GroupData(motifnum).DeficitEditCutoff - Coeff(2) * ced)/Coeff(1);
  	ar = find(ycutoff >= 0);
  	ced = ced(ar);
  	ycutoff = ycutoff(ar);
  	cutofffifty = ycutoff - 0.5*ycutoff(1);
  	ycutoff = min(19.9,ycutoff);
  	cutofffifty = min(19.9,cutofffifty);
  	ar2 = find(cutofffifty >= 0);

    if Grayscale > 0,
			plot(ced,(GroupData(motifnum).DeficitEditCutoff - Coeff(2) * ced)/Coeff(1),'k','linewidth',3);
	    hold on
	    plot(SDR(i,2),SDR(i,1),'k.');
	    plot(SDR(ii,2),SDR(ii,1),'x','color',0.6*[1 1 1]);
	    plot(SDR(iii,2),SDR(iii,1),'x','Color',0.3*[1 1 1],'MarkerSize',8,'LineWidth',2);
	  else
	  	mediumorchid1 = [224	102	255]/255;
	  	thistle1 = [255	225	255]/255;

	  	Lightgreen = [233 250 233]/255;
			Darkgreen  = [211 227 211]/255;

			background = Lightgreen;
			cutofflines = Darkgreen;

			background = thistle1;
			cutofflines = thistle1*0.7;

			patch([0.01 ced 0.01],[0.01 ycutoff 0.01],background)
			hold on
			plot(ced,ycutoff,'color',cutofflines,'linewidth',3);
			plot(ced,cutofffifty,'color',cutofflines,'linewidth',3);
			plot([0 5.5],[0 0],'k');
			plot([0 0],[0 20],'k');
			if ced(ar(end)) < 5.3,
				text(ced(ar(end)),1,'0','color',cutofflines,'fontweight','bold','fontsize',20)
			end
			text(ced(ar2(end)),1,'50','color',cutofflines,'fontweight','bold','fontsize',20)
			text(0.1,1,'100','color',cutofflines,'fontweight','bold','fontsize',20)
	    scatter(SDR(i,2),SDR(i,1),4,'k','filled');
	    scatter(SDR(ii,2),SDR(ii,1),4,'r','filled');
	    plot(SDR(iii,2),SDR(iii,1),'bx','MarkerSize',14,'LineWidth',3);
	  end
		NameForTitle = strrep(GroupData(motifnum).MotifID,'_','\_');
		xlabel(['Minimum interior edit distance to 3D instances in ' NameForTitle],'fontsize',tfs);
		ylabel('Alignment score deficit','fontsize',tfs)
		title(['Random and alignment sequences for ' NameForTitle],'fontsize',tfs);
		set(gca,'fontsize',tfs);
		set(gca,'XTick',0:5);

		axis([0 5.5 0 20]);                           % show every case on the same scale

		if SaveDeficitCoreEditPlot > 0,
			figure(1)
			saveas(gcf,[MSCOutputPath filesep CurrentMotif '_alignment_random_sequences_method_' num2str(GroupData(motifnum).CutoffMethod) '_NumNT_' num2str(GroupData(motifnum).NumNT) '.png']);
			saveas(gcf,[MSCOutputPath filesep CurrentMotif '_alignment_random_sequences_method_' num2str(GroupData(motifnum).CutoffMethod) '_NumNT_' num2str(GroupData(motifnum).NumNT) '.pdf']);
		end

		if PlotCutoffPlane > 0 && MixedScoreCutoff < Inf,
			figure(1)
			hold on
			ax = axis;
			xmax = min(5.5,max(ax(2),2));
			xmax = ax(2);

			xmax = 5.5;

			for x = 0:xmax,
				xx = x * ones(1,1001);                      % core edit values
				yy = ax(3) + (0:1000)*(ax(4)-ax(3))/1000;   % sequence length values
				zz = (MixedScoreCutoff - Coeff(2)*xx)/Coeff(1);        % values that exactly match the cutoff
				ii = find((zz < max(SD(:,1))) .* (zz > min(SD(:,1))));    % values to display on the screen without going above/below the current limits
				ii = find((zz < ax(4)) .* (zz > ax(3)));    % values to display on the screen without going above/below the current limits
				plot3(xx(ii),yy(ii),zz(ii),'r');
			end
			for y = ax(3):ax(4),
				xx = 0 + (0:1000)*ax(2)/1000;               % core edit values
				yy = y * ones(1,1001);                      % sequence length values
				zz = (MixedScoreCutoff - Coeff(2)*xx)/Coeff(1);
				ii = find((zz < max(SD(:,1))) .* (zz > min(SD(:,1))));    % values to display on the screen without going above/below the current limits
				ii = find((zz < ax(4)) .* (zz > ax(3)));    % values to display on the screen without going above/below the current limits
				plot3(xx(ii),yy(ii),zz(ii),'r');
			end
			view([0 0]);
%			axis(ax);
			drawnow
		end

		fprintf('Sensitivity %6.2f%%, Specificity %6.2f%%, Minimum %6.2f%% using method %d\n', 100*TP, 100*TN, 100*min(TP,TN),GroupData(motifnum).CutoffMethod);
		fprintf('Number of false positives with core edit > 0 is %d\n',sum((mynewclasses == 1) .* (Source(k,:) == 2)));
		fprintf('%d * Deficit + %d * Core Edit <= %0.4f\n', Coeff(1), Coeff(2), GroupData(motifnum).DeficitEditCutoff);
		fprintf('Motif index %d\n', iii);
	  fprintf('\n');

	  if pauseafter > 0,
	  	pause
	  end
	end
end

if SaveGroupData > 0,
	save([OutputPath filesep loopType '_GroupData_with_full_cutoffs.mat'],'GroupData');
end

% load([OutputPath filesep loopType '_GroupData_with_full_cutoffs.mat']);

% -------------- Summarize the type of decisions that were made

CM = cat(1,GroupData.CutoffMethod);

for j = 1:max(CM),
	datacounter(j) = sum(CM == j);
end

%datacounter(2) = sum(datacounter([2 4]));
%datacounter(3) = sum(datacounter([3 5]));
%datacounter(6) = sum(datacounter([6 7 8]));

CMText{1} = 'models got the default cutoff from model size';
CMText{2} = 'models had their cutoff set by maximizing TP+TN';
CMText{3} = 'models got the default cutoff plus 2';
CMText{4} = 'models with cutoffs from TP+TN had the cutoff tightened to reduce false positives';
CMText{5} = 'models with default plus 2 had the cutoff tightened to reduce false positives';
CMText{6} = 'models had cutoff set from random sequences only';
CMText{7} = 'random cutoff models had their cutoff made more generous';
CMText{8} = 'random cutoff models had their cutoff made more restrictive';
CMText{9} = 'models had no random sequences and so no model-specific cutoff imposed';
CMText{10} = 'models had the cutoff set from alignment sequences only';
CMText{11} = 'models got the minimum cutoff';

for j = 1:length(CMText),
	fprintf('%3d (%6.2f%%) %s\n', datacounter(j), 100*datacounter(j)/length(GroupData), CMText{j});
end

fprintf('%d groups, total in this table is %d\n',length(GroupData),sum(datacounter));

% -------------- Write out model-specific cutoffs

if isfield(GroupData,'ScoreEditCutoff'),
	clear TN
	clear TP
	for m = 1:length(GroupData),
		GD = GroupData(m);
		TP(m) = GroupData(m).TruePositiveRate;
		TN(m) = GroupData(m).TrueNegativeRate;
		fprintf('Group %3d, %-11s has acceptance rules %8.4f - AlignmentScore <= %8.4f, CoreEdit <= %d, and %8.4f * (%8.4f - AlignmentScore) + %8.4f * CoreEdit <= %8.4f, method %2d,',m, GD.MotifID, max(GD.OwnScore), GD.DeficitCutoff, GD.CoreEditCutoff, GD.DeficitCoeff, max(GD.OwnScore), GD.CoreEditCoeff, GD.DeficitEditCutoff, GD.CutoffMethod);
		fprintf('TP %8.2f%%, TN %8.2f%%, min %8.2f%%, %3d 3D sequences, %5d alignment sequences, %5d random sequences, %4d random matches, %2d NTs, %s\n',100*GD.TruePositiveRate, 100*GD.TrueNegativeRate, 100*min(GD.TruePositiveRate,GD.TrueNegativeRate), GD.NumInstances, GD.NumberOfAlignmentSequences, GD.NumberOfRandomSequences, GD.NumberOfFalsePositives, GD.NumNT, GD.Signature{1});
	end

	if 0 > 1,
		for m = 1:length(GroupData),
			fprintf('Group %3d, %-11s has acceptance rules AlignmentScore >= %8.4f, CoreEdit <= %d, and %8.4f * CoreEdit - %8.4f * AlignmentScore <= %8.4f\n',m, GroupData(m).MotifID, GroupData(m).MinScore, GroupData(m).CoreEditCutoff, GroupData(m).CoreEditCoeff, GroupData(m).DeficitCoeff, GroupData(m).ScoreEditCutoff);
		end
	end

	MSCOutputPath = [OutputPath filesep 'ModelSpecificCutoffs'];
	fid = fopen([MSCOutputPath filesep loopType '_' MotifRelease '_acceptance_rule_numbers.txt'],'w');
	for m = 1:length(GroupData),
		fprintf(fid,'%s\t%0.8f\t%d\t%0.8f\t%0.8f\t%0.8f\n', GroupData(m).MotifID, GroupData(m).MinScore, GroupData(m).CoreEditCutoff, GroupData(m).CoreEditCoeff, GroupData(m).DeficitCoeff, GroupData(m).ScoreEditCutoff);
	end
	fclose(fid);

	for m = 1:length(GroupData),
		fid = fopen([ModelPath filesep GroupData(m).MotifID '_cutoffs.txt'],'w');
		fprintf(fid,'%0.8f\t%d\t%0.8f\t%0.8f\t%0.8f\t%0.8f\t%0.8f\n', GroupData(m).MinScore, GroupData(m).CoreEditCutoff, GroupData(m).CoreEditCoeff, GroupData(m).DeficitCoeff, GroupData(m).ScoreEditCutoff, GroupData(m).DeficitEditCutoff, GroupData(m).MaxScore);
		fclose(fid);
	end

	figure(5)
	clf
	%hist(100*ma(ma > 0),30);
	subplot(2,1,1)
	hist(100*TP(TP > 0),30);
	hist(100*TP,30);
	xlabel(['True positive rate for ' loopType ' models']);
	subplot(2,1,2)
	hist(100*TN(TN > 0),30);
	hist(100*TN,30);
	xlabel(['True negative rate for ' loopType ' models']);

	saveas(gcf,[MSCOutputPath filesep loopType '_TP_TN_rate_histogram.png']);

	figure(3)
	clf
	cutoffmethods = cat(1,GroupData.CutoffMethod);
	scatter(cat(1,GroupData.NumNT),cat(1,GroupData.DeficitEditCutoff),10,cutoffmethods,'filled');
	hold on
	x = 4:18;
	switch loopType,
	case 'IL'
		plot(x,10+1.8*(x-5),'r');
	case 'HL'
		plot(x,10+3*(x-5),'r');
	end
	xlabel('Number of nucleotides');
	ylabel('MixedScore cutoff');
	title('Colored by CutoffMethod');
	caxis([min(cutoffmethods)-0.5 max(cutoffmethods)+0.5])
	colorbar('eastoutside')

	figure(6)
	clf
	scatter(cat(1,GroupData.RandomQuantile96),cat(1,GroupData.DeficitEditCutoff),10,cat(1,GroupData.CutoffMethod),'filled');
	hold on
	plot(0:35,0:35,'r');
	xlabel('96th percentile of MixedScore of random sequences');
	ylabel('MixedScore cutoff');
	title('Colored by CutoffMethod');
	caxis([min(cutoffmethods)-0.5 max(cutoffmethods)+0.5])
	colorbar('eastoutside')

end

diary off
