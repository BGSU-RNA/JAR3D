% pMakeSCFGModels loads all motif groups, makes JAR3D models, writes
% out a fasta file, and writes out all of the correspondences

% MotifLibraryPath is the path to where the .mat files are located
% Release LoopType is like HL/1.8 or IL/1.8; it must match the end of a URL like http://rna.bgsu.edu/rna3dhub/motifs/release/IL/1.6
% OutputPath is the path to where files should be written, to which Release will be appended for organization
% DiagnosticMode is for making different types of models. 1 makes normal models

% Assuming that the Matlab working directory is Motifs, this may work:
% pMakeSCFGModels([pwd filesep 'IL_20140329_0627' filesep 'mat'],pwd,'IL/1.13',1)
% pMakeSCFGModels([pwd filesep '' filesep 'mat'],pwd,'HL/1.13',1)

function [void] = pMakeSCFGModels(MotifLibraryPath,OutputBase,Release,Mode)

if nargin < 4,
  Mode = 1;    % make normal models
end

% ----------------------------------- user controls

MakeEmpiricalDistribution = 0; % no longer making percentile scores

switch Mode,
case 1                      % normal mode to make models.  do not edit
  SampleSize = 20000; % number of samples for empirical distn
  JAR3DSampleSize = 10000;  % maximum number of sequences to pass to JAR3D at a time

  GenerateNewSequences = 0;

  WriteCorrespondences = 1; % use Java to determine correspondences

  WriteModel = 1; % set to 1 to write the models that are made, 0 for testing

  WriteOwnScores = 1;
case 2                      % debugging mode, feel free to edit
  SampleSize = 20000; % number of samples for empirical distn
  JAR3DSampleSize = 10000;  % maximum number of sequences to pass to JAR3D at a time

  GenerateNewSequences = 0;

  WriteCorrespondences = 0; % use Java to determine correspondences

  WriteModel = 0; % set to 1 to write the models that are made, 0 for testing

  WriteOwnScores = 0;
end

JAR3DSampleSize = min(10000,SampleSize);  % maximum number of sequences to pass to JAR3D at a time

loopNum = SampleSize/JAR3DSampleSize;
if loopNum ~= floor(loopNum),
  SampleSize = JAR3DSampleSize * floor(loopNum);
  fprintf('Sample size being adjusted to %d\n',SampleSize);
end

% ----------------------------------- set model type or loop through types

Param = [0 2 0 4 100 1 1 1 100 7 5];                   % See below

% Parameters stored in Param:
% Param(1) verbose
% Param(2) method to use for basepair isostericity
% Param(3) recognize extensible helices and model them as such
% Param(4) adjust substitution probabilities for long-range interactions
% Param(5) how far to look ahead for local basepair interactions
% Param(6) use near interactions
% Param(7) treat insertions as conserved bases
% Param(8) normalize scores for insertions and basepairs:
% Param(9) controls balance between isodiscrepancy and counts; L/(L+Param(9)) is the parameter for counts
% Param(10) tells whether to use BPh interactions and the strength (recommend 7)
% Param(11) tells whether to use BR interactions and the strenghth (recommend 5)

Verbose = Param(1);

Prior = [.5 .5 .5 .5 0];              % Prior distribution for inserted bases

% ------------------------------------------- Set some directories

Release = strrep(Release,'/',filesep);
Release = strrep(Release,'\',filesep);
OutputPath = [OutputBase filesep Release];

loopType = Release(1:2);

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

switch loopType,
case 'JL'
  Rotations = 3;                      % three rotations, for 3-way junctions
case 'IL'
  Rotations = 2;                      % two rotations are computed
case 'HL'
  Rotations = 1;                      % only one "rotation"
end

% --------------------------------- read .mat filenames, select right ones

Filenames = dir(MotifLibraryPath);

keep = [];                               % of all models, which to keep

for m = 1:length(Filenames),
  Filenames(m).modeled = 0;
  if (length(Filenames(m).name) > 2),
    if (Filenames(m).name(1:2) == loopType),
      keep(m) = 1;
      Filenames(m).name = strrep(Filenames(m).name,'.mat','');
    end
  end 
end

Filenames = Filenames(find(keep));

% ----------------------------------- set paths, make directories if needed

DiagnosticBase = [OutputPath filesep 'diagnostic'];
DiagnosticPath = [OutputPath filesep 'diagnostic'];
ModelPath = [OutputPath filesep 'lib'];
ModelDataPath = [OutputPath filesep 'lib' filesep 'data'];
InteractionPath = ModelPath;
SequencePath = ModelPath;

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

if ~(exist(ModelDataPath) == 7) && MakeEmpiricalDistribution > 0,
  mkdir(ModelDataPath);
end

if ~(exist(SequencePath) == 7),
  mkdir(SequencePath);
end

if ~(exist(InteractionPath) == 7),
  mkdir(InteractionPath);
end

Temp.A = 1234;                             % only used for error catching
Temp.B = 9876;                             % only used for error catching

% ----------------------------------- start log file

LogFile = [ModelPath filesep 'log ' date '.txt'];
delete(LogFile);
clc
diary(LogFile);

if MakeEmpiricalDistribution > 0,
  fprintf('Generating %d random sequences for each group\n', SampleSize);
else
  fprintf('Not generating random sequences\n');
end

fprintf('Normalization variable is %d\n', Param(8));

fprintf('pIsoScore2(1,1,4,2) is:\n');

load PairExemplars

pIsoScore2(1,1,4,ExemplarIDI,Param(8))

fprintf('pIsoScore2(1,1,4,2) normalization:\n');

sum(sum(pIsoScore2(1,1,4,ExemplarIDI,Param(8))))

fprintf('Parameter vector is \n');
Param

% ------------------------------------- Gather sequence transition data

TransitionFile = [OutputPath filesep 'transitions.mat'];

if ~(exist(TransitionFile) == 2),
  fprintf('No transitions.mat file found, calculating transition probabilities\n');
  pCalculateTransitions;
else
  fprintf('transitions.mat file found for generating random sequences\n');
end

% ----------------------------------- loop through motifs

clear GroupData                     % structured variable to store data
clear BPhData
clear BRData
BPhcount = 0;                       % number of conserved BPh so far
BRcount  = 0;                       % number of conserved BR so far

clear OwnScoreStats
tic
for m = 1:length(Filenames),
  toc
  tic
  MotifName = Filenames(m).name;

  fprintf('\n');
  disp(['Timestamp ' datestr(now)])

  % ---------- set names of files

  fprintf('pMakeSCFGModels: Analyzing motif group %s\n', MotifName);

  FastaFile = [SequencePath filesep MotifName '.fasta'];
  FastaFileNoGap = [SequencePath filesep MotifName '_nogap.fasta'];
  InteractionFile = [InteractionPath filesep MotifName '_interactions.txt'];
  SubsProbsFile = [ModelPath filesep MotifName '_subsprobs.txt'];
  ModelFile = [ModelPath filesep MotifName '_model.txt'];
  CorrespondenceFile = [ModelPath filesep MotifName '_correspondences.txt'];
  DiagnosticFile = [DiagnosticPath filesep MotifName '_diagnostics.txt'];
  EmpDistFile = [ModelPath filesep MotifName '_distribution.txt'];
  DataFile = [ModelPath filesep MotifName '_data.txt'];
  OwnScoreFile = [ModelPath filesep MotifName '_ownscores.txt'];
  LengthDistFile = [ModelDataPath filesep MotifName '_lengthdist.mat'];
  RandomSequenceDataFile = [ModelDataPath filesep MotifName '_randomsequencedata.mat'];

  AllFileList = [ModelPath filesep 'all.txt'];
  StructuredFileList = [ModelPath filesep 'structured.txt'];
  HTMLMotifList = [ModelPath filesep 'all.html'];

  % ---------- load motif and determine size and number of instances

  clear Search

  load([MotifLibraryPath filesep MotifName '.mat']);

  [L,N] = size(Search.Candidates);
  N = N - 1;

  Filenames(m).numinstances = L;
  Filenames(m).numconserved = N;

  % ---------- temporary fix of base-ribose annotations, until the
  % ---------- RNA 3D motif atlas annotates them correctly

  Search.File = zBaseRiboseInteractions(Search.File);

  % ---------- make an SCFG model of the current model type, write to file

  [Search,Node] = pMakeSingleJAR3DModel(Search,Param,Prior,loopType);

 if isempty(Node),
%    mkdir([MotifLibraryPath filesep 'trouble']);
%    movefile([MotifLibraryPath filesep MotifName '.mat'],[MotifLibraryPath filesep 'trouble' filesep MotifName '.mat']); 
    fprintf('@@@@@@@@@@@@ Motif %s could not be modeled for some reason\n', MotifName);
    Filenames(m).modeled = 0;

 else
    Filenames(m).modeled = 1;

    Text = pNodeToSCFGModelText(Node,5);

    if WriteModel > 0,
      fid = fopen(ModelFile,'w');
      for i = 1:length(Text),
        fprintf(fid,'%s\n', Text{i});
      end
      fclose(fid);
    end

    % ---------- extract sequences of instances into .fasta file

    [Text,T3,T4,T5,ModelFASTA,T7] = xFASTACandidates(Search.File,Search,0,MotifName);

    fid = fopen(FastaFile,'w');
    for t = 1:length(Text),
      fprintf(fid,'%s\n',Text{t});
    end
    fclose(fid);

    % ---------- list out conserved interactions

    Text = pConservedInteractionList(Search,0);

    fid = fopen(InteractionFile,'w');
    for t = 1:length(Text),
      fprintf(fid,'%s\n',Text{t});
    end
    fclose(fid);

    % ---------- list out subs probs in conserved interactions

    Text = pConservedInteractionList(Search,1);

    fid = fopen(SubsProbsFile,'w');
    for t = 1:length(Text),
      fprintf(fid,'%s\n',Text{t});
    end
    fclose(fid);

    % ---------- store interaction signature and other basic data

    GroupData(m).MotifID = MotifName;
    switch loopType,
    case 'HL'    
      GroupData(m).Signature{1} = Search.Signature;
    case 'IL'
      GroupData(m).Signature{1} = Search.Signature;
      GroupData(m).Signature{2} = Search.RSignature;
    end      
    GroupData(m).NumNT = N;
    E = abs(fix(triu(Search.Edge)));
    GroupData(m).NumBasepairs = full(sum(sum((E > 0) .* (E < 14))));

    if GroupData(m).NumBasepairs > Rotations,
      GroupData(m).Structured = 1;            % consider this motif to be "structured"
    else
      GroupData(m).Structured = 0;
    end

    GroupData(m).NumStacks    = full(sum(sum((E > 19) .* (E < 24))));
    BPh = Search.BPh;
    BR  = Search.BR;
    for i = 1:length(BPh(1,:)),
      BPh(i,i) = 0;                      % remove self interactions
      BR(i,i) = 0;
    end
    GroupData(m).NumBPh = full(sum(sum(BPh > 0)));
    GroupData(m).NumBR  = full(sum(sum(BR > 0)));
    GroupData(m).NumInstances = length(Search.Candidates(:,1));
    GroupData(m).Truncate = Search.Truncate;

    si = strrep(GroupData(m).Signature{1},'L','R');
    GroupData(m).NumFixed = length(strfind(si,'-R-'));

    % --------------------------------------- write out base-backbone interactions

    for a = 1:N,
      for b = 1:N,
        if a~= b && Search.BPh(a,b) > 0,
          BPD.CInter = Search.BPh(a,b);             % conserved interaction
          BPD.BCode  = [];                          % base codes
          BPD.Inters = [];                          % individual interactions

          L = length(Search.Candidates(:,1));       % number of candidates
          for c = 1:L,
            f = Search.Candidates(c,N+1);           % file number
            i = Search.Candidates(c,a);
            j = Search.Candidates(c,b);
            BPD.BCode  = [BPD.BCode Search.File(f).NT(i).Code];
            BPD.Inters = [BPD.Inters Search.File(f).BasePhosphate(i,j)];
          end
          BPhcount = BPhcount + 1;
          BPhData(BPhcount) = BPD;

          if L > 1 && Verbose > 0,
            BPD
            for ii = 1:L,
              fprintf('%s ', zBasePhosphateText(BPhData(BPhcount).Inters(ii)));
            end
            fprintf('\n');
          end
        end

        if a~= b && Search.BR(a,b) > 0,
          BRD.CInter = Search.BR(a,b);             % conserved interaction
          BRD.BCode  = [];                          % base codes
          BRD.Inters = [];                          % individual interactions

          L = length(Search.Candidates(:,1));       % number of candidates
          for c = 1:L,
            f = Search.Candidates(c,N+1);           % file number
            i = Search.Candidates(c,a);
            j = Search.Candidates(c,b);
            BRD.BCode  = [BRD.BCode Search.File(f).NT(i).Code];
            BRD.Inters = [BRD.Inters Search.File(f).BaseRibose(i,j)];
          end
          BRcount = BRcount + 1;
          BRData(BRcount) = BRD;

          if L > 1 && Verbose > 0,
            BRD
            for ii = 1:L,
              fprintf('%s ', zBaseRiboseText(BRData(BRcount).Inters(ii)));
            end
            fprintf('\n');
          end
        end
      end
    end

    % ---------- Write motif data file

    fid = fopen(DataFile,'w');
    if strcmp(loopType,'IL'),
      fprintf(fid,'%s %s\n', Search.Signature, Search.RSignature);
    elseif strcmp(loopType,'HL'),
      fprintf(fid,'%s\n', Search.Signature);
    end
    fprintf(fid,'%d nucleotides\n',GroupData(m).NumNT);
    if strcmp(loopType,'IL'),
      fprintf(fid,'%s %s\n', Search.Phonetic, Search.RPhonetic);
    elseif strcmp(loopType,'HL'),
      fprintf(fid,'%s\n', Search.Phonetic);
    end

    fprintf(fid,'%d basepairs\n',GroupData(m).NumBasepairs);
    fprintf(fid,'%d stacks\n',GroupData(m).NumStacks);
    fprintf(fid,'%d base-phosphate\n',GroupData(m).NumBPh);
    fprintf(fid,'%d base-ribose\n',GroupData(m).NumBR);
    fprintf(fid,'%d instances\n',GroupData(m).NumInstances);

    fclose(fid);
    Filenames(m).signature = Search.Signature;

    % ---------- calculate scores of sequences against their own model

    if WriteOwnScores > 0,

      if Verbose > 0,
        fprintf('pMakeSCFGModels: Calculating scores of sequences against their own model\n');
      end

      OwnScores = edu.bgsu.rna.jar3d.JAR3DMatlab.MotifParseSingle(OutputPath,FastaFile,ModelFile);

      % ----------- The following lines prevent the program from being stopped
      % ----------- by a crazy Matlab bug.  It is intermittent, but after a call
      % ----------- to JAR3D, it is hell bent on saying
      % ----------- "Dot name reference on non-scalar structure"

      % Note:  printing Temp.A to the screen is successful when the motif group is a singleton
      % Have a look at OwnScores; there is something about 1 versus many sequences, perhaps

      try
        x = Temp.A + 1;
         Temp.A
      catch ME
        Temp.B = 9876;
      end

      OwnScoreStats{m} = OwnScores;

      fid = fopen(OwnScoreFile,'w');
      for t = 1:length(OwnScores),
        fprintf(fid,'%s\t%16.14f\n',ModelFASTA(t).Sequence,OwnScores(t));
        GroupData(m).OwnScore(t) = OwnScores(t);
        GroupData(m).OwnSequence{t} = ModelFASTA(t).Sequence;
      end
      fclose(fid);
    end

    % ---------- set coefficients for mixed score

    GroupData(m).DeficitCoeff = 1;
    GroupData(m).CoreEditCoeff = 3;

    % ---------- make an empirical distribution for random sequences against model

    if MakeEmpiricalDistribution > 0,

      if exist(RandomSequenceDataFile,'file') && GenerateNewSequences == 0,
        load(RandomSequenceDataFile);
      else
        Scores = zeros(SampleSize,1);
        [Sequences,LengthDist,RandomFASTA] = pMakeRandomSequencesWeighted(Node,loopType,SampleSize,TransitionFile,0);

        save(LengthDistFile,'LengthDist');

        fprintf('pMakeSCFGModels:  Random sequences: ')
        for i = 1:10,
          fprintf('%s  ',RandomFASTA(i).Sequence);
        end
        fprintf('\n');

        loopNum = floor(SampleSize/JAR3DSampleSize);

        fprintf('pMakeSCFGModels:  Parsing random sequences against their own model\n ')

        for i = 1:loopNum,
          SFN = [ModelPath filesep 'RS_' MotifName '_' int2str(i) '.fasta'];
          zWriteFASTA(SFN,RandomFASTA((JAR3DSampleSize*(i-1)+1):(JAR3DSampleSize*i)));
          subScores = edu.bgsu.rna.jar3d.JAR3DMatlab.MotifParseSingle(OutputPath,SFN,ModelFile);

          % ----------- The following lines prevent the program from being stopped
          % ----------- by a crazy Matlab bug.  It is intermittent, but after a call
          % ----------- to JAR3D, it is hell bent on saying
          % ----------- "Dot name reference on non-scalar structure"
          try
            x = Temp.A + 1;
             Temp.A
          catch ME
            Temp.B = 9876;
          end

          Scores(JAR3DSampleSize*(i-1)+1:JAR3DSampleSize*i) = subScores;
          delete(SFN);
        end

        fprintf('Finding core edit distance of %4d randomly-generated sequences against %4d sequences from the motif group\n', length(RandomFASTA), length(ModelFASTA));

        [D1,D2,D3] = pEditDistance(RandomFASTA,ModelFASTA,loopType,'core');

        CoreEditDistance = min(D2,[],2);                      % edit distance in given strand order

        save(RandomSequenceDataFile,'Scores','CoreEditDistance');
      end

%      i = find(CoreEditDistance > 0);                       % sequences with distinct edit distances from 3D instances
%      fprintf('Found and removed %d sequences with core edit distance 0 to a 3D instance, leaving %d sequences\n',length(Scores)-length(i),length(i));
%      CoreEditDistance = CoreEditDistance(i);
%      Scores = Scores(i);

      Deficit = max(GroupData(m).OwnScore) - Scores;
      MixedScores = -(GroupData(m).DeficitCoeff * Deficit + GroupData(m).CoreEditCoeff * CoreEditDistance);

      i = 1:min(2000,length(CoreEditDistance));

      figure(2)
      clf
      plot(CoreEditDistance(i)+0.3*rand(size(CoreEditDistance(i))),Deficit(i),'k.');
      hold on
      plot(0,max(GroupData(m).OwnScore)-GroupData(m).OwnScore,'x','Color',[0.4 0.4 0.4],'MarkerSize',10,'LineWidth',3);   % sequences from 3D
      ax = axis;
      axis([0 ax(2) 0 ax(4)]);
%      axis([0 max(ax([2 4])) 0 max(ax([2 4]))]);
%        axis([0 5 0 20]);

%        plot([0 5],[20 20-5*GroupData(m).CoreEditCoeff],'r');      % line with correct slope

      NameForTitle = strrep(GroupData(m).MotifID,'_','\_');

      xlabel('Minimum core edit distance');
      ylabel('Alignment score deficit');
      title(['Randomly-generated sequences compared to ' NameForTitle]);

      print(gcf,'-dpng',[DiagnosticPath filesep MotifName '_randomsequencescatter.png']);

%        pause


      figure(3)
      clf
      hist(-MixedScores,30);
      h = findobj(gca,'Type','patch');
      set(h,'facecolor',[0.4 0.4 0.4]);
      ax = axis;
      axis([0 ax(2) 0 ax(4)]);

      xlabel(['Deficit + ' num2str(GroupData(m).CoreEditCoeff) ' * CoreEdit']);
      ylabel('Frequency');
      title(['Distribution of randomly-generated sequences for ' NameForTitle]);

      print(gcf,'-dpng',[DiagnosticPath filesep MotifName '_randomsequencehistogram.png']);

      % -------------------------------------- calculate distribution

      Values = unique(MixedScores);           % unique values of mixed score

      N = histc(MixedScores,Values);          % number of occurrences of each unique mixed score

      step = 1/length(MixedScores);           % increase in probability in empirical distribution for each datapoint
      dist = zeros(length(Values),1);
      dist(1) = N(1)*step;
      for i = 2:length(Values)
        dist(i) = dist(i-1) + N(i) * step;
      end

      i = find(dist >= 0.9);
      Precision = 4;                          % number of decimal places in percentile values
      dist(i) = round(dist(i)*(10^Precision))/(10^Precision);
      i = find((dist < 0.9) .* (dist >= 0.8));
      Precision = 3;                          % number of decimal places in percentile values
      dist(i) = round(dist(i)*(10^Precision))/(10^Precision);
      i = find(dist < 0.8);
      Precision = 2;                          % number of decimal places in percentile values
      dist(i) = round(dist(i)*(10^Precision))/(10^Precision);

      prev = -10000;

      fid = fopen(EmpDistFile,'w');
      for i = 1:length(dist)
        cur = dist(i);
        if cur > prev,
          fprintf(fid,'%f %.4f\n',Values(i),dist(i));
        end
        prev = cur;
      end
      fclose(fid);
      clear dist Values MixedScores Sequences

    end

    % ---------- determine correspondences between motif, model, sequences

    T2 = pAlignMotifGroupToModel(Search,Node,MotifName,MotifName); 

    T6 = pColumnsForModel(Node,MotifName);

    clear T8
    for i = 1:length(GroupData(m).OwnScore),
%      T8{i} = sprintf('%s_Instance_%d has_score %0.16f', MotifName, i, GroupData(m).OwnScore(i));
      T8{i} = sprintf('%s_Instance_%d has_score .', MotifName, i);  % more accurate to leave this blank, because the score we have is how the sequence is aligned to the model, which is not right for all sequences!
    end

    T = [T2 T3 T4 T5 T6 T7 T8];

    fid = fopen(CorrespondenceFile,'w');
    for r = 1:length(T),
      T{r} = strrep(T{r},'___','_');
      T{r} = strrep(T{r},'__','_');
      T{r} = strrep(T{r},'has_name _','has_name ');
      fprintf(fid,'%s\n',T{r});
    end
    fclose(fid);

    % ---------- write out correspondences for alignment diagnostics

    if WriteCorrespondences > 0,
      corresp = edu.bgsu.rna.jar3d.JAR3DMatlab.ModelCorrespondences(FastaFile,ModelFile,MotifName,0);
      % ----------- The following lines prevent the program from being stopped
      % ----------- by a crazy Matlab bug.  It is intermittent, but after a call
      % ----------- to JAR3D, it is hell bent on saying
      % ----------- "Dot name reference on non-scalar structure"
      try
        x = Temp.A + 1;
        Temp.A
      catch ME
        Temp.B = 9876;
      end

      corresp = char(corresp);
      correspcell = zStringSplit(corresp,char(10));
      for i = 1:length(correspcell),
        if length(correspcell{i}) > 5,
          correspcell{i} = [GroupData(m).MotifID '_' correspcell{i}];          % prefix with motif ID for this particular diagnostic
          correspcell{i} = strrep(correspcell{i},'___','');
          correspcell{i} = strrep(correspcell{i},'__','');
          correspcell{i} = strrep(correspcell{i},'has_name _','has_name ');
        end
      end

      T = [T correspcell];

      fid = fopen(DiagnosticFile,'w');
      for r = 1:length(T),
        fprintf(fid,'%s\n',T{r});
      end
      fclose(fid);

    end
  end                        % if isempty(Node)
end                           % loop over filenames

% ------------------------------------------------------- write lists of model files

fid = fopen(AllFileList,'w');
for m = 1:length(Filenames),
  if Filenames(m).modeled == 1,
    fprintf(fid,'%s\n',[Filenames(m).name '_model.txt']);
  end
end
fclose(fid);

fid = fopen(StructuredFileList,'w');
for m = 1:length(Filenames),
  if Filenames(m).modeled == 1 && GroupData(m).Structured > 0,
    fprintf(fid,'%s\n',[Filenames(m).name '_model.txt']);
  end
end
fclose(fid);

fprintf('There are %d motif groups, of which %d are structured\n', length(find(cat(1,Filenames.modeled) > 0)), length(find(cat(1,GroupData.Structured)) > 0));

% ------------------------------------------------------- save GroupData

GroupData = GroupData(find(cat(1,Filenames.modeled) > 0));
Filenames = Filenames(find(cat(1,Filenames.modeled) > 0));
save([OutputPath filesep loopType '_GroupData.mat'],'GroupData');

% ------------------------------------------------------- analyze motifs

pMotifCollectionCharacteristics(GroupData);

% ------------------------------------------------------- make HTML list of motifs

fid = fopen(HTMLMotifList,'w');
fprintf(fid,'<html><body><h1>Motifs in release %s</h1>\n',Release);

fprintf(fid,'There are %d motif groups, of which %d are structured<br>\n', length(Filenames), length(find(cat(1,GroupData.Structured)) > 0));

fprintf(fid,'<a href="http://rna.bgsu.edu/rna3dhub/motifs/release/%s">RNA 3D hub site for this release</a>', strrep(Release,filesep,'/'));
fprintf(fid,'<table>\n');
for m = 1:length(Filenames),
  fprintf(fid,'<tr>\n');
  fprintf(fid,'<td><a href="http://rna.bgsu.edu/rna3dhub/motif/view/%s">%s</a></td>\n', Filenames(m).name, Filenames(m).name);

  fprintf(fid,'<td>%d instances</td>\n',Filenames(m).numinstances);
  fprintf(fid,'<td>%d conserved NTs</td>\n',Filenames(m).numconserved);

  fprintf(fid,'<td><a href="%s">Model</a></td>\n',[Filenames(m).name '_model.txt']);
  fprintf(fid,'<td><a href="%s">Interactions</a></td>\n',[Filenames(m).name '_interactions.txt']);
  fprintf(fid,'<td><a href="../diagnostic/%s">Alignment</a></td>\n',[Filenames(m).name '_diagnostics.html']);
  fprintf(fid,'<td>%s</td>\n',Filenames(m).signature);
  fprintf(fid,'<tr>\n');
end
fclose(fid);

diary off

% -------------------------- show statistics about own scores

if Verbose > 0,
  fprintf('pMakeSCFGModels: Statistics on own scores\n');

  for m = 1:length(OwnScoreStats),
    if Verbose > 0,
      sort(OwnScoreStats{m})
      fprintf('\n');
    end

    OwnScoreRange(m) = max(OwnScoreStats{m}) - min(OwnScoreStats{m});
    if OwnScoreRange(m) == 0,
      OwnScoreRange(m) = -5;
    end
  end

  figure(1)
  clf
  hist(OwnScoreRange,30);
end
