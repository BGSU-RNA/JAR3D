% pJAR3DFalsePositiveStudy loops through sequence files corresponding to
% a given motif, runs it against all models, and does diagnostics

% This program accumulates a large amount of data that is useful for exploring the behavior of randomly-generated sequences

function [void] = pJAR3DFalsePositiveStudy(OutputBase,Release,Mode,RandomSequenceMode,numfiles)

if nargin < 4
  RandomSequenceMode = 3;
end

switch Mode
case 1                                    % parse and calculate edit distance, which is slow
  CoreDistSL = 0;                         % minimum sequence length to allow core edit distance matching

  Params.Verbose = 0;

  Params.CutoffType       = 2;            % generic cutoffs
  Params.DeficitCutoff    = 20;           % generic cutoffs
  Params.CoreEditCutoff   = 5;            % generic cutoffs
  Params.PercentileCutoff = 0.2;          % generic cutoffs

  SequenceBySequence = 0;
  CalculateAlignments = 0;                % calculate alignment score, occasionally needed
  CalculateEditDistance = 0;
  PlotPercentileVersusDeficit = 0;
  AccumulateSequenceData = 0;
  AccumulateRandomSequenceData = 0;       % data on sequences for each model
  AccumulateFPData = 0;                   % detailed data about each sequence
  MinimumCutoffScore = -Inf;
  FindFullEditDistance = 0;

  SequenceBySequence = 1;
  Params.CutoffType = 2;                 % use generic cutoffs
  Depth = Inf;                           % accumulate *all* sequences that meet cutoffs, not just the best ones!

case 2                                   % accumulate false positive data
  CoreDistSL = 0;                         % minimum sequence length to allow core edit distance matching

  Params.Verbose = 0;

  Params.CutoffType       = 2;            % generic cutoffs
  Params.DeficitCutoff    = 20;           % generic cutoffs
  Params.CoreEditCutoff   = 5;            % generic cutoffs
  Params.PercentileCutoff = 0.0;          % generic cutoffs; do not impose a percentile cutoff

  SequenceBySequence = 0;
  CalculateAlignments = 0;
  CalculateEditDistance = 0;
  PlotPercentileVersusDeficit = 0;
  AccumulateSequenceData = 0;
  AccumulateRandomSequenceData = 1;       % data on sequences for each model
  AccumulateFPData = 0;                   % detailed data about each sequence
  MinimumCutoffScore = -Inf;
  FindFullEditDistance = 0;

  SequenceBySequence = 1;
  Params.CutoffType = 2;                 % use generic cutoffs
  Depth = Inf;                           % accumulate *all* sequences that meet cutoffs, not just the best ones!

case 3                                    % run the diagnostic using model-specific cutoffs
  CoreDistSL = 0;                         % minimum sequence length to allow core edit distance matching

  Params.Verbose = 1;

  Params.CutoffType       = 3;            % model-specific cutoffs
  Params.DeficitCutoff    = 20;           % generic cutoffs
  Params.CoreEditCutoff   = 5;            % generic cutoffs
  Params.PercentileCutoff = 0.2;          % generic cutoffs

  CalculateAlignments = 0;
  CalculateEditDistance = 0;
  SequenceBySequence = 1;
  PlotPercentileVersusDeficit = 0;
  AccumulateSequenceData = 0;
  AccumulateRandomSequenceData = 0;       % data on sequences for each model
  AccumulateFPData = 0;                   % detailed data about each sequence
  Depth = Inf;
  MinimumCutoffScore = 50;
  FindFullEditDistance = 0;

end

loopType = Release(1:2);

Release = strrep(Release,'/',filesep);
Release = strrep(Release,'\',filesep);
OutputPath = [OutputBase filesep Release];
DiagnosticPath = [OutputBase filesep Release filesep 'diagnostic'];

% create path if it does not already exist
if ~(exist(DiagnosticPath) == 7)
  mkdir(DiagnosticPath);
end

MatchColor(1,:) = 0.6*ones(1,3);
MatchColor(2,:) = 0.3*ones(1,3);
MatchColor(3,:) = 0.0*ones(1,3);

clear SequenceData
clear FPData
clear SD
SDCounter = 0;
FPCounter = 0;

FPModelCounter = zeros(500,1);
SequenceCounter = 0;

ModelTestName = '';                      % will be inserted into paths
DiagnosticTestName = '';

% ----------------------------------- set paths, make directories if needed

if ~(exist(DiagnosticPath) == 7),        % if directory doesn't yet exist
  mkdir(DiagnosticPath);
end

Temp.A = 1234;                         % only used for error catching
Temp.B = 9876;                         % only used for error catching

% ----------------------------------------

fs = 14;                                 % font size for figures
tfs = 13;                                % font size for text
MaxSeqLength = 30;                       % maximum total length to send to JAR3D; maybe not used

if Mode > 2
  figure(1)
  clf
  figure(2)
  clf
end

% --------------------------------- Determine number of rotations

switch loopType
  case 'HL'
    Rotations = 1;                      % only one "rotation"
  case 'IL'
    Rotations = 2;                      % two rotations are computed
  otherwise
    Rotations = str2num(strrep(loopType,'J',''));  % J3, J4, etc.
end

ModelPath = [OutputBase filesep Release filesep 'lib'];
InteractionPath = ModelPath;

% ---------------------------------------- Read data from models

GroupData = pGetModelData(OutputPath,loopType);

if Mode == 2
  for g = 1:length(GroupData)
    GroupData(g).DeficitEditData = [];     % remove any existing data to make room for new
  end
end

if ~isfield(GroupData,'MinScore') && Params.CutoffType == 3
  Params.CutoffType = 2;
  fprintf('Using generic cutoffs because model-specific cutoffs are not yet defined\n');
end

for i = 1:length(GroupData)
  GroupMaxScore(1,i) = max(GroupData(i).OwnScore);
end

Depth = min(Depth,length(GroupData));

% ------------------------------- Tally sequence data from 3D ... seems to not be used

for m = 1:length(GroupData)
  SeqNames{m} = [GroupData(m).MotifID '.fasta'];
end

[FASTA, OwnMotif] = pConcatenateFASTASequences(ModelPath, SeqNames);

clear seqlengths
clear nummodelsbylengths

% switch loopType
% case 'HL'
%   for s = 1:length(FASTA)
%     seqlengths(s,1) = length(FASTA(s).Sequence);
%   end
%   [uniqueseqlengths,i,j] = unique(seqlengths,'rows');

%   maxlengths = max(seqlengths);
%   nummodelsbylengths(maxlengths(1)+3,1) = 0;    % make enough space

%   for k = 1:length(i)
%     a = uniqueseqlengths(k,1);
%     nummodelsbylengths(a,1) = length(unique(OwnMotif(find(j==k))));
%   end

%   newlengths = [];

%   for k = 1:length(uniqueseqlengths(:,1))
%     newlengths(end+1,:) = uniqueseqlengths(k,:) + 1;
%   end

%   [uniqueseqlengths,i,j] = unique([uniqueseqlengths; newlengths],'rows');

%   T = [seqlengths OwnMotif];
%   T = sortrows(T,[1 2]);
% case 'IL'
%   for s = 1:length(FASTA)
%     i = strfind(FASTA(s).Sequence,'*');
%     j = length(FASTA(s).Sequence) - i;
%     seqlengths(s,1) = min(i-1,j);
%     seqlengths(s,2) = max(i-1,j);
%   end
%   [uniqueseqlengths,i,j] = unique(seqlengths,'rows');

%   maxlengths = max(seqlengths);
%   nummodelsbylengths(maxlengths(1)+3,maxlengths(2)+3) = 0;    % make enough space

%   for k = 1:length(i)
%     a = uniqueseqlengths(k,1);
%     b = uniqueseqlengths(k,2);
%     nummodelsbylengths(a,b) = length(unique(OwnMotif(find(j==k))));
%   end

%   newlengths = [];

%   for k = 1:length(uniqueseqlengths(:,1))
%     if uniqueseqlengths(k,1) < uniqueseqlengths(k,2)
%       newlengths(end+1,:) = uniqueseqlengths(k,:) + [1 0];
%     end
%     newlengths(end+1,:) = uniqueseqlengths(k,:) + [0 1];
%     newlengths(end+1,:) = uniqueseqlengths(k,:) + [1 1];
%   end

%   [uniqueseqlengths,i,j] = unique([uniqueseqlengths; newlengths],'rows');

%   T = [seqlengths OwnMotif sum(seqlengths,2)];
%   T = sortrows(T,[1 2 3]);

% % sort by total length???
% %  T = sortrows(T,[4 1 2 3]);
% %  T = T(:,1:3);
% otherwise
%   fprintf('pJAR3DFalsePositiveStudy: loopType %s not recognized\n', loopType);
%   pause
% end

for seqfilenumber = 1:numfiles                           % loop through files of random sequences

  % SVN = [loopType '_FalsePositiveRateTest_' num2str(seqfilenumber) '.fasta'];   % randomly-generated sequences
  SVN = [loopType '_RandomMotifSequences_' num2str(seqfilenumber) '.fasta'];   % randomly-generated sequences

  switch RandomSequenceMode
  case 1
    SVN = [loopType '_RandomMotifSequencesNonCan_' num2str(seqfilenumber) '.fasta'];   % randomly-generated sequences
  case 2
    SVN = [loopType '_RandomMotifSequencesUAOK_' num2str(seqfilenumber) '.fasta'];   % randomly-generated sequences
  case 3
    SVN = [loopType '_RandomMotifSequences_' num2str(seqfilenumber) '.fasta'];   % randomly-generated sequences
  end

  % ---------------------------------------------------------------------

  clear OwnPercentile
  clear SeqLength
  clear LeaveOneOut
  clear NumBetterScore
  clear OwnMLP
  clear OwnTotalProb
  clear SeqNames

  % --------------------------------- Read FASTA file of generated sequences

  CF = strrep(SVN,'.fasta','');
  TCF = strrep(CF,'_','\_');
  CF = strrep(CF,' ','_');

  if Params.Verbose > 0
	  LogFile = [DiagnosticPath filesep CF '_log ' strrep(datestr(now),':','_') '.txt'];
	  diary(LogFile);
	end

  DataFile = [DiagnosticPath filesep CF '.mat'];
  StatusFile = [DiagnosticPath filesep CF '_being_analyzed.txt'];
  Skip = 0;

  if exist(DataFile,'file')
    load(DataFile);
    fprintf('Loaded false positive study data file %s\n',DataFile);
    FirstTime = 0;
  elseif exist(StatusFile,'file') && Mode == 1   % another process is apparently working on this file
    Skip = 1;
  else
    FirstTime = 1;
    fid = fopen(StatusFile,'w');
    fprintf(fid,'Currently being analyzed\n');
    fclose(fid);
    FASTA = zReadFASTA([ModelPath filesep SVN]);
    fprintf('Read random sequence variant file %s\n', SVN);
  end

  % --------------------- Write sequences to one file for each rotation

  for n = 1:length(FASTA)
    FASTA(n).MotifGroup = 'Randomly-generated';
    FASTA(n).Locus = '';
    FASTA(n).Multiplicity = 1;
  end

  NumSequences = length(FASTA);

  [AllSequencesFile,rotatedFASTA] = pWriteSequencesWithRotations(DiagnosticPath,FASTA,loopType,Rotations,CF);

  if Skip == 0 && (CalculateEditDistance > 0 || FirstTime > 0) && length(FASTA) > 0

    % SeqGroup = ones(NumSequences,1) * 0;
    % OwnMotif = ones(NumSequences,1) * 0;

    % ----------------------- Calculate edit distance between FASTA and 3D instances
    CoreEditDistance = Inf * ones(length(FASTA),length(GroupData),length(rotatedFASTA));
    FullEditDistance = Inf * ones(length(FASTA),length(GroupData),length(rotatedFASTA));

    disp(['pJAR3DFalsePositiveStudy Timestamp ' datestr(now)])
    tic;
    fprintf('Finding core edit distance of %4d sequences against %4d models over %4d rotations\n', NumSequences, length(GroupData), length(rotatedFASTA));

    for m = 1:length(GroupData)
      % get all sequences in this motif group
      ModelFASTA = zReadFASTA([ModelPath filesep GroupData(m).MotifID '.fasta']);

      for rotation = 1:length(rotatedFASTA)
        D = pEditDistanceAll(rotatedFASTA{rotation},ModelFASTA,'core');
        CoreEditDistance(:,m,rotation) = min(D,[],2);     % min over sequences in this motif group

        D = pEditDistanceAll(rotatedFASTA{rotation},ModelFASTA,'full');
        FullEditDistance(:,m,rotation) = min(D,[],2);     % min over sequences in this motif group
      end
      if mod(m,40) == 0
        fprintf('pJAR3DFalsePositiveStudy: Checked edit distances for %4d models so far, total time estimate %8.4f minutes\n', m, length(GroupData)*toc/(60*m));
      end
    end

    % if strcmp(loopType,'HL') || strcmp(loopType,'IL')
    %   for m = 1:length(GroupData)
    %     ModelFASTA = zReadFASTA([ModelPath filesep GroupData(m).MotifID '.fasta']);

    %     % D1 is best over all rotations, D2 for original, D3 for first rotation
    %     [D1,D2,D3] = pEditDistance(FASTA,ModelFASTA,loopType,'core');

    %     CoreEditDistance(:,m,1) = min(D2,[],2);     % best match of core seq
    %     if strcmp(loopType,'IL')
    %       CoreEditDistance(:,m,2) = min(D3,[],2);     % best match of core rev'd
    %     end

    %     if mod(m,100) == 0
    %       fprintf('pJAR3DFalsePositiveStudy: Checked interior edit distance for %4d models so far, total time estimate %8.4f minutes\n', m, length(GroupData)*toc/(60*m));

    %     end
    %   end

    %   if FindFullEditDistance > 0
    %     tic;
    %     fprintf('\nFinding full edit distance of %4d sequences against %4d models\n', NumSequences, length(GroupData));
    %     for m = 1:length(GroupData)
    %       ModelFASTA = zReadFASTA([ModelPath filesep GroupData(m).MotifID '.fasta']);

    %       [D1,D2,D3] = pEditDistance(FASTA,ModelFASTA,loopType,'full');

    %       FullEditDistance(:,m,1) = min(D2,[],2);     % best match of full seq
    %       if strcmp(loopType,'IL')
    %         FullEditDistance(:,m,2) = min(D3,[],2);     % best match of full rev'd
    %       end

    %       if mod(m,100) == 0
    %         fprintf('pJAR3DFalsePositiveStudy: Checked full edit distance for %4d models so far, total time estimate %8.4f minutes\n', m, length(GroupData)*toc/(60*m));
    %       end
    %     end
    %   else
    %     FullEditDistance = 99 * ones(size(CoreEditDistance));
    %   end
    % end
  end

	if Skip == 0 && (CalculateAlignments > 0 || FirstTime > 0) && length(FASTA) > 0
    % OwnMotif maps from sequence to MotifName index
    % SeqGroup maps from sequence to 1:L, where L is the number of sequence groups
    %   Each different sequence group gets a different number
    % SeqNames maps from 1:L to a text string for the sequence group

    disp(['pJAR3DFalsePositiveStudy Timestamp ' datestr(now)])

    % ----------------------- Parse all sequences against all models, all rotations

    MLPS       = Inf * ones(length(FASTA),length(GroupData),length(rotatedFASTA));

    % pre-allocating more just caused Matlab to freeze up, don't know why
    % TotalProb  = zeros(length(FASTA),length(GroupData),length(rotatedFASTA));
    % Percentile = zeros(length(FASTA),length(GroupData),length(rotatedFASTA));

    % MLPS(a,m,r) is the score of sequence a against model m using rotation r
    % Percentile(a,m,r) is the percentile of sequence a against model m using r

    tic;
    fprintf('Parsing %d sequences from file %s against %d models over %d rotations\n', NumSequences, SVN, length(GroupData), length(rotatedFASTA));
    for m = 1:length(GroupData)
      for r = 1:Rotations
        % quantileFile = [ModelPath filesep GroupData(m).MotifID '_distribution.txt'];
        MN = [ModelPath filesep GroupData(m).MotifID '_model.txt'];

        S = edu.bgsu.rna.jar3d.JAR3DMatlab.MotifParseSingle(pwd,AllSequencesFile{r},MN);
        MLPS(:,m,r) = S;          % max log probability score for each sequence
%          T = edu.bgsu.rna.jar3d.JAR3DMatlab.MotifTotalProbSingle(pwd,AllSequencesFile{r},MN);
%          MixedScore = -(GroupData(m).DeficitCoeff * (max(GroupData(m).OwnScore) - S) + GroupData(m).CoreEditCoeff * CoreEditDistance(:,m,r));
%          Q = edu.bgsu.rna.jar3d.JAR3DMatlab.getQuantilesFromFile(MixedScore,quantileFile);
%          TotalProb(:,m,r) = T;     % total probability score for each sequence
%          Percentile(:,m,r) = Q;    % percentile of this score
        % TotalProb(:,m,r) = zeros(size(S));     % total probability score for each sequence
        % Percentile(:,m,r) = zeros(size(S));    % percentile of this score
      end
      if mod(m,50) == 0
        fprintf('Parsed against %4d models so far, total time estimate %8.4f minutes\n', m, length(GroupData)*toc/(60*m));
      end
    end
  end

  Marker = 2;
  if Skip == 0 && (FirstTime > 0 || CalculateAlignments > 0 || CalculateEditDistance > 0) && length(FASTA) > 0
    % save(DataFile,'FASTA','Marker','MLPS','TotalProb','Percentile','CoreEditDistance','FullEditDistance');
    save(DataFile,'FASTA','Marker','MLPS','CoreEditDistance','FullEditDistance');

    if exist(StatusFile) > 0
      delete(StatusFile);
    end
  end

  NumSequences = length(FASTA);

  if length(FASTA) > 0
    if Skip == 0 && Mode > 1

      disp(['pJAR3DFalsePositiveStudy Timestamp ' datestr(now)])
      fprintf('Evaluating how sequences fit models\n');

      % --------------------- Evaluate whether each sequence meets the cutoff for each model

      CutoffMet = zeros(size(MLPS));
      CutoffScore = zeros(size(MLPS));
      for mm = 1:length(GroupData)
        for r = 1:Rotations
          Features = [MLPS(:,mm,r) CoreEditDistance(:,mm,r)];
          [CutoffMet(:,mm,r), CutoffScore(:,mm,r)] = pModelSpecificCutoff(GroupData(mm),Features,Params);
        end
      end
    end

    if Skip == 0 && Mode == 2 && AccumulateRandomSequenceData > 0
      MaxScore = zeros(size(MLPS));
      for g = 1:length(GroupData)
        MaxScore(:,g,:) = max(GroupData(g).OwnScore);
      end

      Deficit = max(0,MaxScore - MLPS);          % some might be negative; better score than own model, a strange occurrence

      MixedScore = Deficit + 3*CoreEditDistance;

      for g = 1:length(GroupData)
        [y,br] = min(MixedScore(:,g,:),[],3);      % find the rotation which minimizes the mixed score; this is the best rotation
        for r = 1:Rotations
          i = find((br == r) .* (Deficit(:,g,r) <= Params.DeficitCutoff) .* (CoreEditDistance(:,g,r) <= Params.CoreEditCutoff));
          GroupData(g).DeficitEditData = [GroupData(g).DeficitEditData; [Deficit(i,g,r) CoreEditDistance(i,g,r)]];
        end
      end
    end

    if Skip == 0 && Mode > 2
      % --------------------- Evaluate whether each sequence meets the cutoff for each model

      GenericCutoffMet = zeros(size(MLPS));
      PP = Params;
      PP.CutoffType = 2;                              % generic cutoff, for comparison
      for mm = 1:length(GroupData)
        for r = 1:Rotations
          Features = [MLPS(:,mm,r) CoreEditDistance(:,mm,r)];
          GenericCutoffMet(:,mm,r) = pModelSpecificCutoff(GroupData(mm),Features,PP);
        end
      end

      % --------------------- Determine sequence length(s)

      clear SeqLength
      clear StrandLengths

      SeqLength = zeros(NumSequences,1);
      StrandLengths = zeros(NumSequences,Rotations+1);

      for s = 1:NumSequences
        SeqLength(s) = length(FASTA(s).Sequence)-Rotations+1; % HL, R=1, IL, R=2, etc.

        % StrandLengths is the total sequence length followed by individual strand lengths
        s_parts = strsplit(FASTA(s).Sequence,'*');

        StrandLengths(s,1) = SeqLength(s);
        for k = 1:length(s_parts)
          StrandLengths(s,k+1) = length(s_parts{k});
        end

        % switch loopType
        % case 'IL'
        %   j = strfind(FASTA(s).Sequence,'*');
        %   StrandLengths(s,:) = [length(FASTA(s).Sequence)-1 j-1 length(FASTA(s).Sequence)-j];
        % case 'HL'
        %   StrandLengths(s,:) = [length(FASTA(s).Sequence) length(FASTA(s).Sequence) 1];
        % end
      end

      [y,i] = sort(SeqLength);
      SeqLength = SeqLength(i);
      FASTA = FASTA(i);
      MLPS = MLPS(i,:,:);
      % TotalProb = TotalProb(i,:,:);
      % Percentile = Percentile(i,:,:);
      CoreEditDistance = CoreEditDistance(i,:,:);
      FullEditDistance = FullEditDistance(i,:,:);
      CutoffMet = CutoffMet(i,:,:);
      CutoffScore = CutoffScore(i,:,:);
      GenericCutoffMet = GenericCutoffMet(i,:,:);

      StrandLengths = StrandLengths(i,:);

      MaxScore = zeros(size(MLPS));
      for g = 1:length(GroupData)
      	MaxScore(:,g,:) = max(GroupData(g).OwnScore);
      end

      Deficit = max(0,MaxScore - MLPS);          % some might be negative; better score than own model, a strange occurrence

      % ---------------------------------------- plot percentile versus deficit
      % take the 10 best percentile scores and plot with their deficit
      % and also the 10 best deficit scores and plot with their percentile.

      if PlotPercentileVersusDeficit > 0

        AllPercentiles = [];
        AllDeficits = [];
        AllCoreEdits = [];
        AllSeqLengths = [];

        for s = 1:length(SeqLength)
        	switch loopType
            case 'HL'
             [y,i] = sort(Percentile(s,:,1),2,'descend');                % maximum over rotations, record which rotation
             NewPercentiles = [Percentile(s,i(1:Depth),1)];
             NewDeficits    = [Deficit(s,i(1:Depth),1)];
             NewEdits       = [CoreEditDistance(s,i(1:Depth),1)];
             [y,i] = sort(NewPercentiles,2,'descend');

             AllPercentiles = [AllPercentiles NewPercentiles(i(1:Depth))];
             AllDeficits    = [AllDeficits NewDeficits(i(1:Depth))];
             AllCoreEdits   = [AllCoreEdits NewEdits(i(1:Depth))];
             AllSeqLengths  = [AllSeqLengths ones(1,Depth)*SeqLength(s)];

             [y,i] = sort(Deficit(s,:,1),2,'ascend');                % maximum over rotations, record which rotation
             NewPercentiles = [Percentile(s,i(1:Depth),1)];
             NewDeficits    = [Deficit(s,i(1:Depth),1)];
             NewEdits       = [CoreEditDistance(s,i(1:Depth),1)];
             [y,i] = sort(NewDeficits,2,'ascend');

             AllPercentiles = [AllPercentiles NewPercentiles(i(1:Depth))];
             AllDeficits    = [AllDeficits NewDeficits(i(1:Depth))];
             AllCoreEdits   = [AllCoreEdits NewEdits(i(1:Depth))];
             AllSeqLengths  = [AllSeqLengths ones(1,Depth)*SeqLength(s)];
           case 'IL'
    			 	[y,i] = sort(Percentile(s,:,1),2,'descend');                % maximum over rotations, record which rotation
    	    	[z,j] = sort(Percentile(s,:,2),2,'descend');                % maximum over rotations, record which rotation
    			 	[y,i] = sort(MLPS(s,:,1),2,'descend');                % maximum over rotations, record which rotation
    	    	[z,j] = sort(MLPS(s,:,2),2,'descend');                % maximum over rotations, record which rotation
    	    	NewMLPS        = [MLPS(s,i(1:Depth),1) MLPS(s,j(1:Depth),2)];
    	    	NewPercentiles = [Percentile(s,i(1:Depth),1) Percentile(s,j(1:Depth),2)];
    	    	NewDeficits    = [Deficit(s,i(1:Depth),1)    Deficit(s,j(1:Depth),2)];
    	    	NewEdits       = [CoreEditDistance(s,i(1:Depth),1)  CoreEditDistance(s,j(1:Depth),2)];
    	    	NewFullEdits   = [FullEditDistance(s,i(1:Depth),1)  FullEditDistance(s,j(1:Depth),2)];
    	    	ModelNumbers   = [i(1:Depth) j(1:Depth)];
    	    	[y,k] = sort(NewPercentiles,2,'descend');
    	    	[y,k] = sort(NewMLPS,2,'descend');

    	    	AllPercentiles = [AllPercentiles NewPercentiles(k(1:Depth))];
    	    	AllDeficits    = [AllDeficits NewDeficits(k(1:Depth))];
    	    	AllCoreEdits   = [AllCoreEdits NewEdits(k(1:Depth))];
    	    	AllSeqLengths  = [AllSeqLengths ones(1,Depth)*SeqLength(s)];

            k = k(1:Depth);

            k = find((Percentile(s,:,1) > 0.5) + (Percentile(s,:,2) > 0.5) + (Deficit(s,:,1) < 50) + (Deficit(s,:,2) < 50) + (CoreEditDistance(s,:,1) <= 10) + (CoreEditDistance(s,:,2) <= 10));

            if AccumulateSequenceData > 0
    		    	for d = 1:length(k)
    		        SD.Percentile   = NewPercentiles(k(d));
    		        SD.Deficit      = NewDeficits(k(d));
    		        SD.CoreEdit     = NewEdits(k(d));
    		        SD.FullEdit     = NewFullEdits(k(d));
    		        SD.Conserved    = GroupData(ModelNumbers(k(d))).NumNT;
    		        SD.NumInstances = GroupData(ModelNumbers(k(d))).NumInstances;
    		        SD.MLPS         = NewMLPS(k(d));
    		        SD.Length       = SeqLength(s);
    		        SD.Source       = 2 + (d-1)*0.2;
    		        SD.Source       = 2;
    	%	        SD.Sequence     = FASTA(s).Sequence;
    	%	        SD.LoopID       = '';
    		        SD.MotifID      = GroupData(ModelNumbers(k(1))).MotifID;

    		        if mod(SDCounter,1000) == 0
    		          SequenceData(SDCounter + 1000) = SD;
    		        end

    		        SDCounter = SDCounter + 1;
    		        SequenceData(SDCounter) = SD;
    		      end
    		    end

    	    	if 0 > 1
    		    	[y,i] = sort(Deficit(s,:,1),2,'ascend');                % maximum over rotations, record which rotation
    		    	[z,j] = sort(Deficit(s,:,2),2,'ascend');                % maximum over rotations, record which rotation
    		    	NewPercentiles = [Percentile(s,i(1:Depth),1) Percentile(s,j(1:Depth),2)];
    		    	NewDeficits    = [Deficit(s,i(1:Depth),1)    Deficit(s,j(1:Depth),2)];
    		    	NewEdits       = [CoreEditDistance(s,i(1:Depth),1)  CoreEditDistance(s,j(1:Depth),2)];
    		    	[y,i] = sort(NewDeficits,2,'ascend');

    		    	AllPercentiles = [AllPercentiles NewPercentiles(i(1:Depth))];
    		    	AllDeficits    = [AllDeficits NewDeficits(i(1:Depth))];
    		    	AllCoreEdits   = [AllCoreEdits NewEdits(i(1:Depth))];
    		    	AllSeqLengths  = [AllSeqLengths ones(1,Depth)*SeqLength(s)];
    		    end
          otherwise
            fprintf('pJAR3DFalsePositiveStudy: loopType %s not recognized\n', loopType);
    	    end
        end

        if 0 > 1
    	    for j = unique(SeqLength)
    		    k = find(AllSeqLengths == j);
    	%	    k = find((AllSeqLengths == j) .* (AllCoreEdits >= 2));

    		    figure(j)
    		    scatter(AllDeficits(k),100*AllPercentiles(k),4,min(6,AllCoreEdits(k)),'filled');
    		    hold on
    		    axis([0 15 0 100]);
    		    caxis([0 7])
    				xlabel('Alignment score deficit; colored by core edit distance, max 6');
    				ylabel('Percentile against models');
    				title(['Best ' num2str(Depth) ' scores for each random ' loopType ' sequence of length ' num2str(j)]);
    				if seqfilenumber == 2
    					colorbar('eastoutside')
    				end
    			end
    		end
      end

  		if 0 > 1
  		 	figure(10);
  		 	clf
  	    for g = 1:length(GroupData)
  				scatter(Deficit(:,g,1),100*Percentile(:,g,1),4,min(FullEditDistance(:,g,1),1),'filled');
  				hold on
  				scatter(Deficit(:,g,2),100*Percentile(:,g,2),4,min(FullEditDistance(:,g,2),1),'filled');
  				xlabel('Alignment score deficit');
  				ylabel('Percentile against own model');
  				title(['Scores of random sequences against model with ' num2str(GroupData(g).NumNT) ' nucleotides']);
  				axis([0 15 0 100.001]);
  			end
  		end

  		% --------------------------------------- sequence by sequence analysis

  		if SequenceBySequence > 0
  	    [b,i,j] = unique(StrandLengths,'rows');
  	    % [y,k] = sortrows(b,[1 2]);
  	    % [y,k] = sortrows(b,[2 3]);                        % sort by shorter strand length, then longer
        [s,t] = size(b);
        [y,k] = sortrows([b sum(b,2)],[t+1 1]);           % sort by decreasing total length

        LengthsToPosition = containers.Map();       % code each set of strand lengths by a number
        Match = zeros(length(k),2,2);
        FullMatch = zeros(length(k),2);
        FullEditOne = zeros(length(k),2);
        CoreEditOne = zeros(length(k),2);
        CoreEditZero = zeros(length(k),2);
        CoreMatch = zeros(length(k),2);
        GoodMatch = zeros(length(k),2);
        TotalSequences = zeros(length(k),2);

  	    for w = 1:length(k)
          LengthsToPosition(mat2str(y(w,1:(Rotations+1)))) = w;
          % fprintf('%s\n',mat2str(y(w,1:(Rotations+1))));
  	    end
  	    LengthPairs = y;

  	    for n = 1:length(FASTA)
  	      SL = SeqLength(n);                                             % plain sequence length
          % fprintf('%s\n',mat2str(StrandLengths(n,:)));
  	      sl = LengthsToPosition(mat2str(StrandLengths(n,:)));  % map strand lengths to position

          % StrandLengths(n,:)
          % fprintf('sequence %s position %d\n', FASTA(n).Sequence,sl);

  	      TotalSequences(sl,1) = TotalSequences(sl,1) + 1;    % count total number of sequences with these lengths

  	      [m,r] = min(CoreEditDistance(n,:,:),[],3);          % minimum over rotations and record which rotation
  	      [MinEditCore,mm] = min(m);                          % minimum over models and record which model
  	      r = r(mm);                                          % optimal rotation at optimal model

  	      edindex = MinEditCore + 1;                          % to store the minimum edit distance

          if edindex > size(TotalSequences,1)
            % awkward but works
            Match(edindex,2) = 0;
            FullMatch(edindex,2) = 0;
            FullEditOne(edindex,2) = 0;
            CoreEditOne(edindex,2) = 0;
            CoreEditZero(edindex,2) = 0;
            CoreMatch(edindex,2) = 0;
            GoodMatch(edindex,2) = 0;
            TotalSequences(edindex,2) = 0;
          end

  	      TotalSequences(edindex,2) = TotalSequences(edindex,2) + 1;    % count total number of sequences

  	      if MinEditCore == 0                                 % exact core sequence match
  	        CoreMatch(sl,1) = CoreMatch(sl,1) + 1;
  	        CoreMatch(edindex,2) = CoreMatch(edindex,2) + 1;
  	        cmmm = mm;
  	        cmr = r;
  	      end

          if MinEditCore == 0                                 % close to Core sequence match
            CoreEditZero(sl,1) = CoreEditZero(sl,1) + 1;
            CoreEditZero(edindex,2) = CoreEditZero(edindex,2) + 1;
          end

  	      if MinEditCore <= 1                                 % close to Core sequence match
  	        CoreEditOne(sl,1) = CoreEditOne(sl,1) + 1;
  	        CoreEditOne(edindex,2) = CoreEditOne(edindex,2) + 1;
  	      end

  	      [m,r] = min(FullEditDistance(n,:,:),[],3);          % minimum over rotations and record which rotation
  	      [MinEditFull,mm] = min(m);                          % minimum over models and record which model
  	      r = r(mm);                                          % optimal rotation at optimal model

  	      if MinEditFull == 0                                 % exact full sequence match
  	        FullMatch(sl,1) = FullMatch(sl,1) + 1;
  	        FullMatch(edindex,2) = FullMatch(edindex,2) + 1;
  	        fmmm = mm;
  	        fmr = r;
  	      end

  	      if MinEditFull <= 1                                 % close to full sequence match
  	        FullEditOne(sl,1) = FullEditOne(sl,1) + 1;
  	        FullEditOne(edindex,2) = FullEditOne(edindex,2) + 1;
  	      end

  	      SequenceCounter = SequenceCounter + 1;

          mm = [];
          rm = [];

          if MinEditFull == 0
          	mm = fmmm;
          	rm = fmr;
          elseif MinEditCore == 0 && SL >= CoreDistSL
          	mm = cmmm;
          	rm = cmr;
          end

          if max(max(GenericCutoffMet(n,:,:))) > 0
            Match(sl,1,2)      = Match(sl,1,2) + 1;                        % generic match
            Match(edindex,2,2) = Match(edindex,2,2) + 1;
          end

          gmm = mm;                                                        % good matches
          grm = rm;

        	[yy,ii] = sort(max(MLPS(n,:,:),[],3),2,'descend');               % best alignment scores first
          [yy,ii] = sort(max(CutoffScore(n,:,:),[],3),2,'descend');        % best cutoff scores first
          jj = 1;
          while length(mm) < Depth && jj <= length(MLPS(n,:,1)),           % loop through all models
          	cmn = ii(jj);                                                  % current model number
            [cut,r] = max(CutoffMet(n,cmn,:));                             % find rotation that meets the cutoff

            if cut > 0
              mm = [mm cmn];
              [cut,r] = max(CutoffScore(n,cmn,:));                         % find rotation that maximizes the score
              rm = [rm r];
            end

            if cut > MinimumCutoffScore
              gmm = [gmm cmn];
              [cut,r] = max(CutoffScore(n,cmn,:));                         % find rotation that maximizes the score
              grm = [grm r];
            end

            jj = jj + 1;
          end

          if length(gmm) > 0
            GoodMatch(sl,1) = GoodMatch(sl,1) + 1;
            GoodMatch(edindex,2) = GoodMatch(edindex,2) + 1;
          end

          if length(mm) > 0
            Match(sl,1,1) = Match(sl,1,1) + 1;                    % match with current parameters
            Match(edindex,2,1) = Match(edindex,2,1) + 1;

            allmm = mm;
            allrm = rm;

            for kk = 1:length(allmm)
            	mm = allmm(kk);
            	rm = allrm(kk);

  		        if CoreEditDistance(n,mm,rm) > 0
  			        FPModelCounter(mm) = FPModelCounter(mm) + 1;
  			      end

  		        if AccumulateFPData > 0
  			        SD.Percentile   = Percentile(n,mm,rm);
  			        SD.Deficit      = max(GroupData(mm).OwnScore)-MLPS(n,mm,rm);
  			        SD.CoreEdit     = CoreEditDistance(n,mm,rm);
  			        SD.FullEdit     = FullEditDistance(n,mm,rm);
  			        SD.Conserved    = GroupData(mm).NumNT;
  			        SD.NumInstances = GroupData(mm).NumInstances;
  			        SD.MLPS         = MLPS(n,mm,rm);
  			        SD.Length       = length(FASTA(n).Sequence)-Rotations+1;
  			        SD.Source       = 2;
  			        SD.Basepairs    = GroupData(mm).NumBasepairs;
                SD.NumFixed     = GroupData(mm).NumFixed;
  			        SD.Sequence     = FASTA(n).Sequence;
  			        SD.MotifID      = GroupData(mm).MotifID;

  			        if mod(FPCounter,1000) == 0
  			          SD(FPCounter + 1000) = SD;
  			        end

  			        FPCounter = FPCounter + 1;
  			        FPData(FPCounter) = SD;
  				    end
  				  end

            mm = allmm(1);
            rm = allrm(1);

            if Params.Verbose > 0
  		        fprintf('%s', FASTA(n).Sequence);
              fprintf(' %d ', CutoffMet(n,mm,rm));
  	          fprintf(' Matches %2d; best is ',length(allmm));
  	          fprintf('%11s, %2d NTs, ', GroupData(mm).MotifID, GroupData(mm).NumNT);
              fprintf(' %20s,',GroupData(mm).OwnSequence{1});
  	          fprintf('score %10.2f, ', MLPS(n,mm,rm));
  	          fprintf('deficit %5.2f, ', max(GroupData(mm).OwnScore)-MLPS(n,mm,rm));
  	          % fprintf('prct %6.2f%%, ', 100*Percentile(n,mm,rm));
              fprintf('CutoffScore %6.2f, ', CutoffScore(n,mm,rm));
  	          fprintf('ed %2d,%2d ', FullEditDistance(n,mm,rm), CoreEditDistance(n,mm,rm));
  	          fprintf('%-40s', GroupData(mm).Signature{1});
    				  fprintf('\n');
  	        end

          else
          	if Params.Verbose > 0
  %		          fprintf(' Has no good match                                              ');
  	        end
          end

  	      if Params.Verbose > 0
  %		      fprintf('\n');
  		    end
  	    end

  	    for v = 1:2
  	      figure(v)
  	      clf
  	      FullMatchperc    = 100 * FullMatch(:,v) ./ TotalSequences(:,v);
  	      FullEditOneperc  = 100 * FullEditOne(:,v) ./ TotalSequences(:,v);
  	      CoreMatchperc    = 100 * CoreMatch(:,v) ./ TotalSequences(:,v);
  	      CoreEditOneperc  = 100 * CoreEditOne(:,v) ./ TotalSequences(:,v);
          CoreEditZeroperc = 100 * CoreEditZero(:,v) ./ TotalSequences(:,v);
          GoodMatchperc    = 100 * GoodMatch(:,v) ./ TotalSequences(:,v);
          Matchperc        = 100 * Match(:,v,1) ./ TotalSequences(:,v);
          GenericMatchperc = 100 * Match(:,v,2) ./ TotalSequences(:,v);

  	      switch v
  	      case 1
  	        range = 1:length(LengthPairs(:,1));
  	        index = 1:length(LengthPairs(:,1));

            range = 1:min(80,length(LengthPairs(:,1)));
            index = range;

  	        subplot(2,1,1);
  	      case 2
  	        range = 0:min(6,length(LengthPairs(:,1)));
  	        index = range + 1;
  	      end

          hold on
  %        plot(range,FullMatchperc(index),'k:','linewidth',2);
  %        plot(range,CoreEditOneperc(index),'k--','linewidth',2);

  %        plot(range,GenericMatchperc(index),'linewidth',2,'color',MatchColor(1,:));

          plot(range,CoreEditZeroperc(index),'linewidth',2,'color',0.7*[1 1 1]);
          plot(range,GoodMatchperc(index),'linewidth',2,'color',0.35*[1 1 1]);
          plot(range,Matchperc(index),'k','linewidth',2);

  %        plot(range,Matchperc(index),'k.','markersize',16);

  %          plot(range,FullMatchperc(index),'r--','linewidth',2);
  %          plot(range,FullMatchperc(index),'r','linewidth',2);
  %            plot(range,FullMatchperc(index),'r.','markersize',16);
  %            plot(range,FullEditOneperc(index),'g--','linewidth',2);
  %            plot(range,FullEditOneperc(index),'g.','markersize',16);
  %            plot(range,CoreMatchperc(index),':','linewidth',2,'color',[148  0 211]/255);
  %            plot(range,CoreMatchperc(index),'.','markersize',16,'color',[148  0 211]/255);
  %          plot(range,CoreEditOneperc(index),'b','linewidth',2);
  %            plot(range,CoreEditOneperc(index),'b.','markersize',16);
  %            plot(range,FullMatchperc(index),'r.','markersize',16);
       %       plot(range,FullEditOneperc(index),'g--','linewidth',2);
       %       plot(range,FullEditOneperc(index),'g.','markersize',16);
  %            plot(range,CoreMatchperc(index),':','linewidth',2,'color',[148  0 211]/255);
  %            plot(range,CoreMatchperc(index),'.','markersize',16,'color',[148  0 211]/255);
  %            plot(range,CoreEditOneperc(index),'b.','markersize',16);

  	      ylabel('Percentage','fontsize',tfs)

  	      switch v
  	      case 1
            title(['Acceptance rate of sequences in ' loopType '\_Rand by strand length'],'fontsize',tfs)
  	        % for w = 1:max(range)
            %   switch loopType
            %   case 'HL'
    	      %     text(w,-1,[sprintf('%d (%d)',nummodelsbylengths(LengthPairs(w,2)), LengthPairs(w,2))],'rotation',90,'horizontalalignment','right','FontName','FixedWidth','fontsize',6);
            %   case 'IL'
            %     text(w,-1,[num2str(nummodelsbylengths(LengthPairs(w,2),LengthPairs(w,3))) '(' num2str(LengthPairs(w,2)) ',' sprintf('%2d',LengthPairs(w,3)) ')'],'rotation',90,'horizontalalignment','right','FontName','FixedWidth','fontsize',6);
            %   end
  	        % end
  	        axis([0 max(range) 0 100.5])
  	        set(gca,'xtick',[])
  	%        xlabel('Sequence lengths')
  	        print(gcf,'-dpng',[DiagnosticPath filesep CF '_RandomSequenceMatchRate_sequence_length.png']);
  	        print(gcf,'-dpdf',[DiagnosticPath filesep CF '_RandomSequenceMatchRate_sequence_length.pdf']);
  	      case 2
            title(['Acceptance rate of sequences in ' loopType '\_Rand by interior edit distance'],'fontsize',tfs)
  	        for w = index
  	          text(range(w),10,num2str(TotalSequences(w,2)),'rotation',90,'horizontalalignment','right');
  	        end
  	        axis([0 25 0 100.1])
  	        xlabel('Minimum interior edit distance to known 3D instances')
  	        print(gcf,'-dpng',[DiagnosticPath filesep CF '_RandomSequenceMatchRate_edit_distance.png']);
  	        print(gcf,'-dpdf',[DiagnosticPath filesep CF '_RandomSequenceMatchRate_edit_distance.pdf']);
  	      end

          % loop over strand lengths
          k = keys(LengthsToPosition);
          fprintf('[total_length strand_lengths] Matchperc GoodMatchperc CoreEditZeroperc');
          for i = 1:length(k)
            sl = LengthsToPosition(k{i});
            if sl > 0 && v == 1
              fprintf('%s\t%8.4f\t%8.4f\t%8.4f\n',k{i},Matchperc(sl),GoodMatchperc(sl),CoreEditZeroperc(sl));
              % fprintf('%d\t%8.4f\t%8.4f\t%8.4f\t%d\n',a,Matchperc(sl),GoodMatchperc(sl),CoreEditZeroperc(sl),nummodelsbylengths(a,1));
            end
          end

          % switch loopType
          % case 'IL'
          %   for a = 1:length(LengthsToPosition(:,1))
          %     for b = a:length(LengthsToPosition(1,:))
          %       sl = LengthsToPosition(a,b);
          %       if sl > 0 && v == 1
          %         fprintf('%d\t%d\t%8.4f\t%8.4f\t%8.4f\n',a,b,Matchperc(sl),GoodMatchperc(sl),CoreEditZeroperc(sl));
          %         % fprintf('%d\t%d\t%8.4f\t%8.4f\t%8.4f\t%d\n',a,b,Matchperc(sl),GoodMatchperc(sl),CoreEditZeroperc(sl),nummodelsbylengths(a,b));
          %         CEZRate(a,b) = CoreEditZeroperc(sl);
          %         GMRate(a,b) = GoodMatchperc(sl);
          %         MRate(a,b) = Matchperc(sl);
          %       end
          %     end
          %   end
          % case 'HL'
          %   % loop over strand lengths
          %   k = keys(LengthsToPosition);
          %   for i = 1:length(k)
          %     sl = LengthsToPosition(k{i});
          %     if sl > 0 && v == 1
          %       fprintf('%s\t%8.4f\t%8.4f\t%8.4f\n',k{i},Matchperc(sl),GoodMatchperc(sl),CoreEditZeroperc(sl));
          %       % fprintf('%d\t%8.4f\t%8.4f\t%8.4f\t%d\n',a,Matchperc(sl),GoodMatchperc(sl),CoreEditZeroperc(sl),nummodelsbylengths(a,1));
          %     end
          %   end
          % end

          % grayscale heat map for IL

          if 0>1 && strcmp(loopType,'IL')
            figure(v+2)
            clf

            colormap(gray)
            map = colormap;
            map = 1-map;

            if v == 1
              subplot(2,2,1)
              pcolor(CEZRate);
              shading flat
              colormap(map);
              colorbar('eastoutside')
              title('Percentage with interior edit distance 0')
              subplot(2,2,2)
              pcolor(GMRate);
              shading flat
              colormap(map);
              colorbar('eastoutside')
              title('Percentage with cutoff score 50 or over')
              subplot(2,2,3)
              pcolor(MRate);
              shading flat
              colormap(map);
              colorbar('eastoutside')
              title('Percentage accepted by at least one motif group')
              subplot(2,2,4)
              [s,t] = size(MRate);
              % pcolor(nummodelsbylengths(1:s,1:t));
              shading flat
              colormap(map);
              colorbar('eastoutside')
              title('Number of motif groups with instance of the given size')
              print(gcf,'-dpng',[DiagnosticPath filesep loopType '_FP_heatmap.png']);
            end
          end
  	    end

  	    fprintf('\n\n\n\n\n');
      end
	  end
	end

  if AccumulateFPData > 0
    FPData = FPData(1:FPCounter);
  % [DiagnosticPath filesep loopType '_Random_Sequence_Data.mat'],'SequenceData');
    save([OutputPath filesep loopType '_Random_Sequence_Data_FP_Only.mat'],'FPData');
    fprintf('Saved false positive data\n');
  end

%  FPModelIDs = cat(1,FPData.MotifID);

  if Mode > 2
    FPModelCounter = FPModelCounter(1:length(GroupData));
    [y,i] = sort(FPModelCounter,1,'ascend');
    for k = 1:length(i)
      GD = GroupData(i(k));
      if isfield(GroupData,'TruePositiveRate')
    		fprintf('Model %-11s matched a sequence %3d times, TP %6.2f%%, TN %6.2f%%, %3d instances, signature %s, sequence %s\n', GD.MotifID,FPModelCounter(i(k)),100*GD.TruePositiveRate,100*GD.TrueNegativeRate,GD.NumInstances,GD.Signature{1},GD.OwnSequence{1});
      else
        fprintf('Model %-11s matched a sequence %3d times, TP        , TN        , %3d instances, signature %s, sequence %s\n', GD.MotifID,FPModelCounter(i(k)),GD.NumInstances,GD.Signature{1},GD.OwnSequence{1});
    	end
    end
    fprintf('Overall false positive rate is %0.4f%%\n', 100*sum(Match(:,1,1))/SequenceCounter);
    fprintf('Number of models giving false positives is %d out of %d\n',sum(FPModelCounter > 0),length(GroupData));
    fprintf('Rate at which sequences have interior edit distance 0 is %8.2f\n',100*sum(CoreEditZero(:,1))/SequenceCounter);
    fprintf('Rate at which sequences have cutoff score greater than 50 is %8.2f\n',100*sum(GoodMatch(:,1))/SequenceCounter);

    diary off

    figure(3)
    clf
    hist(FPModelCounter,30);
    xlabel('Number of false positives');
    ylabel('Number of models');
  end
end

if Mode == 2 && AccumulateRandomSequenceData > 0
  save([OutputPath filesep loopType '_GroupData.mat'],'GroupData');
  fprintf('Saved group data\n');
  GroupData
end

if 0 > 1
	for j = unique(SeqLength)
		figure(j);
		print(gcf,'-dpng',[DiagnosticPath filesep loopType '_Percentile_Deficit_Top_' num2str(Depth) '_Scores_' num2str(j) '_NT.png']);
		print(gcf,'-dpdf',[DiagnosticPath filesep loopType '_Percentile_Deficit_Top_' num2str(Depth) '_Scores_' num2str(j) '_NT.pdf']);
	end
end

% save([DiagnosticPath filesep loopType '_Random_Sequence_Data_Above_Cutoffs.mat'],'SequenceData');
