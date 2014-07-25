% pJAR3DDiagnostics runs JAR3D on all motif sequences against all motif models.

% Release is like HL/1.8 or IL/1.8; it must match the end of a URL like http://rna.bgsu.edu/rna3dhub/motifs/release/IL/1.8
% OutputPath is the path to where files should be written, to which Release will be appended for organization
% SequenceSource is relative to OutputPath.  It is often OutputPath / Release, but could be different
% DiagnosticMode is for making different types of models. 1 makes normal models

function [void] = pJAR3DDiagnostics(OutputBase,Release,SequenceSource,DiagnosticMode)

if nargin < 4,
  DiagnosticMode = 1;
end

PercentileCutoffs = [0.8 0.9 0.95];
DeficitCutoffs = [10 8 4];

Params.SizeOfGuessSet = 1;
Params.UseMultiplicity = 0;
Params.CoreDistSL = 4;
Params.MinMotifSize = 6;
Params.MinSeqLengthIL = 6;
Params.MaxSeqLength = 30;
Params.NumSequencesToShow = Inf;

% basic cutoffs only
Params.CutoffType       = 2;      % use generic cutoffs
Params.DeficitCutoff    = 20;
Params.CoreEditCutoff   = 5;
Params.PercentileCutoff = 0.2;

Verbose = 1;

ModelTestName = '';                      % will be inserted into paths
DiagnosticTestName = '';

% Specify the type of diagnostic with DiagnosticMode
%   example:  DiagnosticMode = 1;       % do nothing special
%   name it with DiagnosticTestName
%   example:  DiagnosticTestName = 'UU at break'; % will be inserted into paths

% ---------------------------------------- set directories

Release = strrep(Release,'\',filesep);
Release = strrep(Release,'/',filesep);
OutputPath = [OutputBase filesep Release];

SequenceSource = strrep(SequenceSource,'\',filesep);
SequenceSource = strrep(SequenceSource,'/',filesep);

if strcmp(Release,SequenceSource),
  UsingSequencesFromMotifGroups = 1;
  SequencePath = [OutputBase filesep SequenceSource filesep 'lib'];
else
  UsingSequencesFromMotifGroups = 0;
  SequencePath = [OutputBase filesep SequenceSource];
end

% ------------------------------------ check for basic user information

DM{1} = '';            % DiagnosticMode 1 use raw sequences
DM{2} = '_UU_at_break';% DiagnosticMode 2 insert UU pair at *
DM{3} = '_UU_at_end';  % DiagnosticMode 3 insert UU pair at beginning/end
DM{4} = '_UU_at_both'; % DiagnosticMode 4 insert UU pairs at * and beginning/end
DM{5} = '_rescore_insertions'; % DiagnosticMode 5 remove insertions and re-score

if ~exist('DiagnosticMode'),
  disp('Please specify a diagnostic mode, for example DiagnosticMode = 1;')
  disp('Using DiagnosticMode = 1');
  DiagnosticMode = 1;
end

% --------------------------------- Determine number of rotations

loopType = Release(1:2);

switch loopType,
case 'JL'
  Rotations = 3;                      % three rotations, for 3-way junctions
case 'IL'
  Rotations = 2;                      % two rotations are computed
case 'HL'
  Rotations = 1;                      % only one "rotation"
end

% ----------------------------------- set paths, make directories if needed

DiagnosticBase = [OutputPath filesep 'diagnostic'];
DiagnosticPath = [OutputPath filesep 'diagnostic'];
ModelPath = [OutputPath filesep 'lib'];
InteractionPath = ModelPath;

if ~(exist(OutputPath) == 7),        % if directory doesn't yet exist
  mkdir(OutputPath);
end

if ~(exist(DiagnosticBase) == 7),
  mkdir(DiagnosticBase);             
end

if ~(exist(DiagnosticPath) == 7),
  mkdir(DiagnosticPath);
end

if ~(exist(ModelPath) == 7),
  mkdir(ModelPath);
end

if ~(exist(SequencePath) == 7),
  mkdir(SequencePath);
end

if ~(exist(InteractionPath) == 7),
  mkdir(InteractionPath);
end

% ---------------------------------------- 

fs = 14;                                 % font size for figures
tfs = 13;

% ---------------------------------------- Read data from models

[GroupData,MotifEquivalence] = pGetModelData(OutputPath,loopType);

% ---------------------------------------------------------------------

clear OwnPercentile
clear SeqLength
clear LeaveOneOut
clear NumBetterScore
clear OwnMLP
clear OwnTotalProb

% --------------------------------- Read sequence file names
% --------------------------------- Store sequence files in FASTA

clear SeqNames

if UsingSequencesFromMotifGroups == 1,     % internal diagnostic
  for m = 1:length(GroupData),
    SeqNames{m} = [GroupData(m).MotifID '.fasta'];
  end
  [FASTA, OwnMotif] = pConcatenateFASTASequences(SequencePath, SeqNames);
  SeqGroup = OwnMotif;                      % same concept here
  for n = 1:length(FASTA),
    FASTA(n).Multiplicity = 1;
  end
else
%  SequencePath = [pwd filesep 'lib' filesep SequenceSource];
  DiagnosticBase = [pwd filesep Release filesep 'externaldiagnostic' filesep SequenceSource DiagnosticTestName];
  DiagnosticPath = [pwd filesep Release filesep 'externaldiagnostic' filesep SequenceSource DiagnosticTestName];
  if ~(exist(DiagnosticPath) == 7),
    mkdir(DiagnosticPath);
  end

%  [FASTA, OwnMotif, SeqGroup, SeqNames] = pReadMSASequences(SequencePath,MotifNames);
	fprintf('pJAR3DDiagnostics:  This option is no longer available\n');

end

% -------------------------------- Calculate core edit distances between instances within each motif group

max1 = 0;
for g = 1:length(GroupData),
  i = find(OwnMotif == g);

  max2 = 0;
  for a = 1:length(i),
    for b = (a+1):length(i),
      d = pEditDistance(FASTA(i(a)),FASTA(i(b)),loopType,'core');
      max1 = max(max1,d);
      max2 = max(max2,d);
    end
  end

  if 0 > 1,                    % look at core edit distances between 3D instances
    fprintf('Motif group %s has maximum core edit distance %2d\n', GroupData(g).MotifID, max2);
    if max2 > 5,
      for j = i,
        fprintf('  %s\n',FASTA(i).Sequence);
      end
    end
  end
end

fprintf('Overall maximum intragroup core edit distance is %2d\n', max1);

% OwnMotif maps from sequence to MotifName index
% SeqGroup maps from sequence to 1:L, where L is the number of sequence groups 
%   Each different sequence group gets a different number
% SeqNames maps from 1:L to a text string for the sequence group

NumSequences = length(FASTA);
NumSeqGroups = max(SeqGroup);

% --------------------- Add UU pairs to FASTA for diagnostic purposes

if any(DiagnosticMode == [2 3 4]),
  FASTA = pAddPairstoFASTAforDiagnostic(FASTA,NumSequences,DiagnosticMode);
end

% --------------------- Write sequences to one file for each rotation

[AllSequencesFile] = pWriteSequencesWithRotations(ModelPath,FASTA,loopType,Rotations);

% ------------------------- Load parsing and edit distance data if necessary

DataFile = [DiagnosticPath filesep 'ParsingData.mat'];

if exist(DataFile,'file'),
	load(DataFile)
	fprintf('Loaded data file %s\n', DataFile);
else
  % ----------------------- Calculate edit distance between FASTA and 3D instances
  clear CoreEditDistance
  clear AvgCoreEditDistance

  fprintf('Finding core edit distance of %4d sequences against %4d models\n', NumSequences, length(GroupData));
  for m = 1:length(GroupData),
    ModelFASTA = zReadFASTA([ModelPath filesep GroupData(m).MotifID '.fasta']); % read sequences from *models*

    [D1,D2,D3] = pEditDistance(FASTA,ModelFASTA,loopType,'core');
    if strcmp(Release,SequenceSource),                % internal diagnostic
%     i = find(OwnMotif == m);                        % seqs from current motif
%     D2(i,:) = zMakeZeroInf(D2(i,:));                % ignore exact matches
%     D3(i,:) = zMakeZeroInf(D3(i,:));                % ignore exact matches
    end
    CoreEditDistance(:,m,1) = min(D2,[],2);     % best match of core seq
    AvgCoreEditDistance(:,m,1) = mean(D2,2);     % best match of core seq
    if strcmp(loopType,'IL'),
      CoreEditDistance(:,m,2) = min(D3,[],2);     % best match of core rev'd
      AvgCoreEditDistance(:,m,2) = mean(D3,2);     % best match of core rev'd
    end

    if mod(m,50) == 0,
      fprintf('Checked core edit distance for %4d models so far\n', m);
    end
  end

  % ----------------------- Calculate edit distance between FASTA and 3D instances
  clear FullEditDistance
  clear AvgFullEditDistance

  fprintf('Finding full edit distance of %4d sequences against %4d models\n', length(FASTA), length(GroupData));
  for m = 1:length(GroupData),
    ModelFASTA = zReadFASTA([ModelPath filesep GroupData(m).MotifID '.fasta']);

    [D1,D2,D3] = pEditDistance(FASTA,ModelFASTA,loopType,'full');
    if strcmp(Release,SequenceSource),           % internal diagnostic
%      i = find(OwnMotif == m);                        % seqs from current motif
%      D2(i,:) = zMakeZeroInf(D2(i,:));                % ignore exact matches
%      D3(i,:) = zMakeZeroInf(D3(i,:));                % ignore exact matches
    end
    FullEditDistance(:,m,1) = min(D2,[],2);     % best match of Full seq
    AvgFullEditDistance(:,m,1) = mean(D2,2);     % best match of Full seq
    if strcmp(loopType,'IL'),
      FullEditDistance(:,m,2) = min(D3,[],2);     % best match of Full rev'd
      AvgFullEditDistance(:,m,2) = mean(D3,2);     % best match of Full rev'd
    end

    if mod(m,50) == 0,
      fprintf('Checked full edit distance for %4d models so far\n', m);
    end
  end

	% ----------------------- Parse all sequences against all models, all rotations

	clear MLPS
	clear TotalProb
	clear Percentile

	% MLPS(a,m,r) is the score of sequence a against model m using rotation r
	% Percentile(a,m,r) is the percentile of sequence a against model m using r

	Temp.A = 12;
	Temp.B = 27;

	fprintf('Parsing %4d sequences against %4d models\n', NumSequences, length(GroupData));
	for m = 1:length(GroupData),
	  quantileFile = [ModelPath filesep GroupData(m).MotifID '_distribution.txt'];
	  for r = 1:Rotations,
	    MN = [ModelPath filesep GroupData(m).MotifID '_model.txt'];
	    S = edu.bgsu.rna.jar3d.JAR3DMatlab.MotifParseSingle(pwd,AllSequencesFile{r},MN);
%	    T = edu.bgsu.rna.jar3d.JAR3DMatlab.MotifTotalProbSingle(pwd,AllSequencesFile{r},MN);

      MixedScore = -(GroupData(m).DeficitCoeff * (max(GroupData(m).OwnScore) - S) + GroupData(m).CoreEditCoeff * CoreEditDistance(:,m,r));

%	    Q = edu.bgsu.rna.jar3d.JAR3DMatlab.getQuantilesFromFile(MixedScore,quantileFile);

	    MLPS(:,m,r) = S;          % max log probability score for each sequence
%	    TotalProb(:,m,r) = T;     % total probability score for each sequence
%	    Percentile(:,m,r) = Q;    % percentile of this score
      TotalProb(:,m,r) = zeros(size(S));     % total probability score for each sequence
      Percentile(:,m,r) = zeros(size(S));    % percentile of this score
	  end

	  if mod(m,50) == 0,
	    fprintf('Parsed against %4d models so far\n', m);
	  end
	end

	% ----------- The following lines prevent the program from being stopped
	% ----------- by a crazy Matlab bug.  It is intermittent, but after a call
	% ----------- to JAR3D, it is hell bent on saying
	% ----------- "Dot name reference on non-scalar structure"
	try
	  Temp.A
	catch ME
	  Temp.B = 27;
	end

  save(DataFile,'FASTA','MLPS','TotalProb','Percentile','CoreEditDistance','FullEditDistance','AvgCoreEditDistance','AvgFullEditDistance');

end

% --------------------- Evaluate whether each sequence meets the cutoff for each model

CutoffMet   = zeros(size(MLPS));
CutoffScore = zeros(size(MLPS));
Par = Params;
Par.CutoffType = 3;
for mm = 1:length(GroupData),
  for r = 1:Rotations,
    Features = [MLPS(:,mm,r) CoreEditDistance(:,mm,r) Percentile(:,mm,r)];
    [CutoffMet(:,mm,r) CutoffScore(:,mm,r)] = pModelSpecificCutoff(GroupData(mm),Features,Params);
  end
end

% --------------------- Determine sequence length, excess length, percentile

clear GroupSize

for g = 1:max(SeqGroup),
  j = find(SeqGroup == g);
  GroupSize(g) = length(j);           % number of sequences in this group
end

clear ExcessSeqLength
clear SeqLength

for s = 1:NumSequences,
  SeqLength(s) = length(FASTA(s).Sequence)-Rotations+1; % HL, R=1, IL, R=2, etc.

  g = OwnMotif(s);
  OwnPercentile(s) = Percentile(s,g,1);
  OwnMLP(s) = MLPS(s,g,1);
  OwnTotalProb(s) = TotalProb(s,g,1);
  OwnEditDistance(s) = CoreEditDistance(s,g,1);
  OwnCutoffScore(s) = CutoffScore(s,g,1);
end

for s = 1:NumSequences,
  j = find(OwnMotif == OwnMotif(s));     % sequences from the same group
  ExcessSeqLength(s) = SeqLength(s) - mode(SeqLength(j));
  OwnDeficit(s) = max(OwnMLP(j)) - OwnMLP(s); 
end


% --------------------- Histogram Cutoff Score against own model

figure(1)
clf
hist(OwnCutoffScore,30)
h = findobj(gca,'Type','patch');
set(h,'facecolor',[0.4 0.4 0.4]);
xlabel('Cutoff Score','fontsize',tfs);
ylabel('Frequency','fontsize',tfs);
title(['Cutoff Score of ' num2str(length(FASTA)) ' ' loopType ' sequences from 3D structures'],'fontsize',tfs);
ax = axis;
set(gca,'fontsize',tfs);
q = 100*[0.99 0.95 0.9 0.8 0.5 0];
for m = 1:length(q),
  text(ax(1)+0.1*(ax(2)-ax(1)), (0.98-0.08*m)*ax(4), sprintf('%5.2f%% score higher than %3.0f', 100*sum(OwnCutoffScore > q(m))/NumSequences, q(m)),'fontsize',tfs);
end
print(gcf,'-dpng',[DiagnosticPath filesep 'Own_CutoffScore_Histogram' DM{DiagnosticMode} '.png']);

% ------------------------------- Group-group diagnostic, all groups
% ------------------------------- Count models with better MLPS and Percentiles

OnlyStructured = 0;

NBS = pGroupGroupDiagnostic(NumSeqGroups,GroupData,OnlyStructured,FASTA,MLPS,OwnMotif,SeqGroup);

figure(1)
clf
T = ['MLP group-group diagnostic, ' num2str(NumSeqGroups) ' sequence groups against ' num2str(length(GroupData)) ' models'];
pHistogramNumBetterScore(NBS,T,fs);

print(gcf,'-dpng',[DiagnosticPath filesep 'Group_Group' DM{DiagnosticMode} '.png']);

% ------------------------------- Structured Group-group diagnostic
% ------------------------------- Count models with better MLPS and Percentiles

OnlyStructured = 1;

NBS = pGroupGroupDiagnostic(NumSeqGroups,GroupData,OnlyStructured,FASTA,MLPS,OwnMotif,SeqGroup);

figure(1)
clf

NumStructured = sum(cat(1,GroupData.Structured)==1);
T = ['MLP structured group-group diagnostic, ' num2str(NumStructured) ' sequence groups against ' num2str(NumStructured) ' models'];
pHistogramNumBetterScore(NBS,T,fs);

print(gcf,'-dpng',[DiagnosticPath filesep 'Group_Group_Structured' DM{DiagnosticMode} '.png']);

% ------------------------------- Individual-group diagnostic

OnlyStructured = 0;

NBS = pIndividualGroupDiagnostic(GroupData,OnlyStructured,OwnMotif,MLPS);

figure(1)
T = ['MLP individual-group diagnostic, ' num2str(NumSequences) ' sequences against ' num2str(length(GroupData)) ' models'];

pHistogramNumBetterScoreIndividual(NBS,abs(ExcessSeqLength),T,fs);
ylabel('Colored by difference in sequence length','fontsize',fs);

print(gcf,'-dpng',[DiagnosticPath filesep 'Individual_Group' DM{DiagnosticMode} '.png']);


figure(2)
clf
T = ['MLP individual-group diagnostic, ' num2str(NumSequences) ' sequences against ' num2str(length(GroupData)) ' models'];
cc = ceil(1 + 10*(cat(1,GroupData.NumBasepairs)-2)./cat(1,GroupData.NumNT));
cc = cc(OwnMotif);
pHistogramNumBetterScoreIndividual(NBS,cc,T,fs);
ylabel('Colored by core basepairs per nucleotide','fontsize',fs)

% --------------------------------------- Individual sequence details
% This can be run after the individual-group or structured individual-group diagnostic

LogFile = [DiagnosticPath filesep 'log ' date '.txt'];
delete(LogFile);
clc
diary(LogFile);

Text = pIndividualGroupSequenceRundown(Params,OnlyStructured,OwnMotif,GroupData,MLPS,FASTA,ModelPath,SeqGroup,OwnEditDistance,CoreEditDistance,Percentile,3,Verbose,CutoffScore,FullEditDistance,AvgCoreEditDistance,CutoffMet,MotifEquivalence);

% write Text to text file

SequenceRundownFile = [DiagnosticPath filesep 'rundown ' strrep(Release,filesep,' ') ' ' date '.txt'];
fid = fopen(SequenceRundownFile,'w');
for i = 1:length(Text),
  fprintf(fid,'%s\n',Text{i});
end
fclose(fid);

diary off

% ------------------------------- Use multiple sequences from each group

pMultipleGroupDiagnostic

% --------------------- Histogram number of instances in each sequence group

clear GroupSize

for g = 1:NumSeqGroups,
  GroupSize(g) = length(find(SeqGroup == g));
end

figure(1)
clf
hist(GroupSize,30);
title('Number of instances in each group','fontsize',fs);
set(gca,'fontsize',fs)
print(gcf,'-dpng',[DiagnosticPath filesep 'Num_Instances' DM{DiagnosticMode} '.png']);

figure(1)
clf
semilogy(1:NumSeqGroups,GroupSize,'.');
title('Number of instances versus group number','fontsize',fs);
set(gca,'fontsize',fs)
print(gcf,'-dpng',[DiagnosticPath filesep 'Num_Instances_By_Group' DM{DiagnosticMode} '.png']);

% ------------------------------- Basic statistics for the paper

NumStructured = sum(cat(1,GroupData.Structured)==1);

NumSingleton = sum(cat(1,GroupData.NumInstances) == 1);

fprintf('%6d %s motif groups\n',length(GroupData),loopType);
fprintf('%6d loop instances\n',length(FASTA));
fprintf('%6d singleton groups\n',NumSingleton);
fprintf('%6d with basepairs or BR or BPh in addition to the flanking cWW pair(s)\n',NumStructured);
