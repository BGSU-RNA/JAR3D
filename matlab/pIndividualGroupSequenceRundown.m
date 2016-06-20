% pIndividualGroupSequenceRundown analyzes the alignment of individual sequences to probabilistic models and reports on the models that the sequence fits best

% Copy Text into a spreadsheet, then sort by columns B, S, and U to get a nice ordering

function [Text,ROC,Confusion,Correctness] = pIndividualGroupSequenceRundown(Params,OnlyStructured,OwnMotif,GroupData,MLPS,FASTA,ModelPath,SeqGroup,OwnCoreEditDistance,CoreEditDistance,Percentile,Criterion,Verbose,CutoffScore,FullEditDistance,AvgCoreEditDistance,CutoffMet,MotifEquivalence)

CoreDistSL = Params.CoreDistSL;
SizeOfGuessSet = Params.SizeOfGuessSet;
UseMultiplicity = Params.UseMultiplicity;

Confusion = zeros(length(GroupData),length(GroupData)+1);   % to count mis-classifications

Correctness.Criterion = Criterion;
Correctness.NumSeqs   = length(FASTA);
Correctness.TotalMultiplicity = sum(cat(1,FASTA.Multiplicity));
Correctness.AcceptedByMultiplicity = 0;
Correctness.AcceptedByMultiplicityExclInteriorMatch = 0;
Correctness.AcceptedByMultiplicityExclFullMatch = 0;
Correctness.CorrectByMultiplicity = 0;
Correctness.CorrectByMultiplicityExclFullMatch = 0;
Correctness.CorrectByMultiplicityExclInteriorMatch = 0;
Correctness.TotalMultiplicityExclFullMatch = 0;
Correctness.TotalMultiplicityExclInteriorMatch = 0;
Correctness.CorrectSequences = 0;
Correctness.NumNT = GroupData(OwnMotif(1)).NumNT;
Correctness.MeanSequenceLength = GroupData(OwnMotif(1)).MeanSequenceLength;
Correctness.GroupNum = OwnMotif(1);
Correctness.MotifID = GroupData(OwnMotif(1)).MotifID;
Correctness.NumBasepairs = GroupData(OwnMotif(1)).NumBasepairs;
Correctness.NumBPh = GroupData(OwnMotif(1)).NumBPh;
Correctness.MostCommonSequence = FASTA(1).Sequence;

if nargin < 17,
  Criterion = 3;
end

Score = 10000*(FullEditDistance == 0);  % full sequence matches always win

switch Criterion,
case 1,                                  % alignment score
  Score = Score + MLPS;
case 2,                                  % edit distance, ties broken by alignment score
  Score = Score - CoreEditDistance + MLPS/1000;
case 3,                                  % plain edit distance, choose randomly
  Score = Score - CoreEditDistance;
case 4,                                  % average edit distance to known instances
  Score = Score + MLPS - CoreEditDistance;
case 5,                                  % total probability
  Score = Score + CutoffScore;             % temporarily not really total probabilty
case 6,                                  % percentile score
  Score = Score + Percentile + MLPS/1000000;    % use MLPS to break ties but not affect the decision more than that
case 7,                                  % just check to see if the cutoff is met, don't really score by Score
  Score = Score + MLPS;
end

PointColors = 'gcbkr';
loopType = GroupData(1).MotifID(1:2);

if nargin < 18,
  Verbose = 1;
end

ROC = [];
ROC(2,2) = 0;                              % avoid plotting problems

NumBetterModelsToShow = 5;

if Verbose == 1,
  if OnlyStructured == 1,
    fprintf('Structured individual-group diagnostic\n');
  else
    fprintf('Individual-group diagnostic, all groups\n');
  end
  fprintf('Listing groups by interaction signature\n');
end

%fprintf('Listing all sequences in all groups for which at least one sequence scores better against another model or has percentile worse than 0.2 against its own model\n');

if OnlyStructured == 1,
  str = find(cat(1,GroupData.Structured) == 1);  % structured internal loops
else
  str = 1:length(GroupData);
end

seqnum = 1;                                         % row of cell array for text output, one line for each sequence
Text{seqnum} = ['1' char(9) '0' char(9) 'Motif ID' char(9) 'NTs - number of nucleotides in this motif group' char(9) 'Instances' char(9) 'BP - number of conserved basepairs, including the flanking cWW pairs' char(9) 'BPh' char(9) 'BR' char(9) 'Structured? - 1 if it has more than two basepairs, 0 if not.  The next column is to record the human annotation.' char(9) 'Is this motif internally structured?' char(9) 'Loop id' char(9) 'Sequence' char(9) 'URL of the motif' char(9) 'Signature' char(9) 'Score - max log probability score' char(9) 'Percentile' char(9) 'Edit distance - core edit distance to a distinct sequence within this motif group' char(9) 'Number of models that scored better' char(9) 'Explanation for there being better scoring models - one word so we can classify using this word!' char(9) 'Comment' char(9) 'ID - of the better-scoring motif' char(9) 'URL' char(9) 'Signature' char(9) 'Score' char(9) 'Percentile' char(9) 'Edit distance' char(9) 'NTs - number of nucleotides in this motif' char(9) 'URL' char(9) 'Signature' char(9) 'Score' char(9) 'Percentile' char(9) 'Edit distance' char(9) 'NTs - number of nucleotides in this motif'];
seqnum = seqnum + 1;

for g = 1:length(GroupData),
	Sig{g} = GroupData(g).Signature{1};
end

[y,grouporder] = sort(Sig);

for gg = 1:length(grouporder),                 % run through sequence groups
  Text{seqnum} = [num2str(seqnum) char(9) num2str(2*gg-1) char(9)];              % blank line between sequence groups
  seqnum = seqnum + 1;

  g = grouporder(gg);

  s = find(SeqGroup == g);                     % sequences in this group

  if length(s) > 0,
    mn = OwnMotif(s(1));                                % correct model number

    alignmentnumrows = 0;
    alignmentzeroedit = 0;

    for k = 1:length(s),
      n = s(k);
      if OnlyStructured == 0 || GroupData(mn).Structured > 0,
        alignmentnumrows = alignmentnumrows + FASTA(n).Multiplicity;

        if OwnCoreEditDistance(n) == 0,
          alignmentzeroedit = alignmentzeroedit + FASTA(n).Multiplicity;
        end
      end
    end

    Correctness.PercentZeroEdit = alignmentzeroedit/alignmentnumrows;

    em = find(MotifEquivalence(mn,:) > 0);     % this model and all equivalent to it

    if Verbose == 1,
      fprintf('Group %3d is from %s which you can see at %s\n', g, GroupData(mn).MotifID, ['http://rna.bgsu.edu/rna3dhub/motif/view/' GroupData(mn).MotifID]);

      if GroupData(mn).Structured == 1,
        fprintf('This group is considered to be structured ***************************\n');
      else
        fprintf('This group is considered to be unstructured ---------------------------\n');
      end

      fprintf('Number of NTs: %2d  Signature: %s\n', GroupData(mn).NumNT, GroupData(mn).Signature{1});
      fprintf('%0.2f%% of rows in the alignment have core edit distance 0 from a sequence in 3D\n', 100*alignmentzeroedit/alignmentnumrows);
      fprintf('Sequence (multiplicity in alignment) percentile, core edit distance\n');

      for ww = 1:length(em),
        if em(ww) ~= mn,
          fprintf('Equivalent motif: %s %s\n',GroupData(em(ww)).MotifID, GroupData(em(ww)).Signature{1});
        end
      end
    end

    % ------------------------------- order sequence matches approximately best to worst

    if length(s) > 1,
      clear OrderMatrix
      for k = 1:length(s),
        n = s(k);                                         % sequence number
        [m,r] = max(Score(n,:,:),[],3);                   % maximum over rotations
        bettermodel = str(find(m(str) >= Score(n,mn,1))); % model numbers of better models
        bettermodel = setdiff(bettermodel,mn);            % take out own model

        OrderMatrix(k,1) = length(bettermodel);
        OrderMatrix(k,2) = round(10000*Percentile(n,mn,1));
        OrderMatrix(k,3) = CoreEditDistance(n,mn,r(mn));
        OrderMatrix(k,4) = MLPS(n,mn,1);
      end
      [y,i] = sortrows(OrderMatrix,[3 -2 -4]);             % sort by edit distance, number with better score, percentile
      s = s(i);
    end

    bestmodels = [];                                      % keep track for evaluation of voting as a classification technique

    for k = 1:min(length(s),Params.NumSequencesToShow),   % loop through sequences
      n = s(k);                                           % sequence number

      if UseMultiplicity > 0,
        Counter = FASTA(n).Multiplicity;                  % use multiplicity of this variant
      else
        Counter = 1;                                      % simply count sequences
      end

      ed = CoreEditDistance(n,mn,1);                      % minimum edit distance to instances known from 3D
      ed = min(ed,30);                                    % for 3D structures, ed can be Inf
      [p,q] = size(ROC);
      if p <= ed,
        ROC(ed+1,2) = 0;                                  % make space for a larger edit distance
      end

      LocalScore = Score(n,:,:);

      SL = length(strrep(FASTA(n).Sequence,'*',''));

      if SL >= CoreDistSL,
        LocalScore = LocalScore + 1000*(CoreEditDistance(n,:,:) == 0);  % core sequence matches
      end

      [m,r] = max(LocalScore,[],3);                       % maximum score over rotations

      OwnCutoffMet = max(CutoffMet(n,em,1));              % 1 if any cutoffs are met

      if OwnCutoffMet > CutoffMet(n,mn,1) && Verbose == 1,
        fprintf('Equivalent motif meets cutoff but not original one\n');
      end

      if OwnCutoffMet == 1,
        OwnScore = max(m(em));       % maximum score over equivalent models
      else
        OwnScore = -Inf;
      end

      if OwnScore > m(mn) && Verbose == 1,
        fprintf('Equivalent motif has better score\n');
      end

      bettermodel = [];                         % keep track of better models
      equalmodel = [];

      for w = 1:length(str),                    % loop through desired models
        cmn = str(w);                           % current model number
        if cmn ~= mn,
          if CutoffMet(n,cmn,r(cmn)) > 0,
            if m(cmn) > OwnScore,                 % better score
              bettermodel = [bettermodel cmn];
            end
            if m(cmn) == OwnScore && MotifEquivalence(mn,cmn) <= 0, % different model scores the same
              equalmodel = [equalmodel cmn];
            end
          end
        end
      end

      numbetter = length(bettermodel);
      numequal  = length(equalmodel);

      if numbetter > 0,
        [y,i] = sort(-m(bettermodel));                    % sort best to worst
        bettermodel = bettermodel(i);                     % re-order
      end

      if SizeOfGuessSet == 1,
        if OwnCutoffMet > 0 && numbetter == 0,
          if numequal > 0,
            Confusion(mn,mn) = Confusion(mn,mn) + 1/(1+numequal);  % tie
            Confusion(mn,equalmodel) = Confusion(mn,numequal) + 1/(1+numequal);
          else
            Confusion(mn,mn) = Confusion(mn,mn) + 1;      % correct match
          end
        elseif numbetter > 0,
          j = find(m(bettermodel) == m(bettermodel(1)));                              % tied better matches
          Confusion(mn,bettermodel(j)) = Confusion(mn,bettermodel(j)) + 1/length(j);  % each one gets a bump
        else
          Confusion(mn,length(GroupData)+1) = Confusion(mn,length(GroupData)+1) + 1; % no match
        end
      end

      score = 0;

      if OwnCutoffMet > 0,
        Correctness.AcceptedByMultiplicity = Correctness.AcceptedByMultiplicity + FASTA(n).Multiplicity;
        if OwnCoreEditDistance(n) > 0,
          Correctness.AcceptedByMultiplicityExclInteriorMatch = Correctness.AcceptedByMultiplicityExclInteriorMatch + FASTA(n).Multiplicity;
        end
        if FullEditDistance(n,mn,1) > 0,
          Correctness.AcceptedByMultiplicityExclFullMatch = Correctness.AcceptedByMultiplicityExclFullMatch + FASTA(n).Multiplicity;
        end
        if Criterion == 7,
          score = 1;                                % does the sequence fit the correct model well enough?
        elseif numbetter + numequal < SizeOfGuessSet,
          score = 1;
        elseif numbetter < SizeOfGuessSet,
          score = min(1,(SizeOfGuessSet - numbetter) / (numequal + 1));
        end
      end

      if Params.CoreEditZeroSuccess > 0 && OwnCoreEditDistance(n) == 0,    % simply recognize these as being correct
        OwnCutoffMet = 1;
        score = 1;
      end

      Correctness.CorrectByMultiplicity = Correctness.CorrectByMultiplicity + FASTA(n).Multiplicity * score;
      if OwnCoreEditDistance(n) > 0,
        Correctness.TotalMultiplicityExclInteriorMatch = Correctness.TotalMultiplicityExclInteriorMatch + FASTA(n).Multiplicity;
        Correctness.CorrectByMultiplicityExclInteriorMatch = Correctness.CorrectByMultiplicityExclInteriorMatch + FASTA(n).Multiplicity * score;
      end
      if FullEditDistance(n,mn,1) > 0,
        Correctness.TotalMultiplicityExclFullMatch = Correctness.TotalMultiplicityExclFullMatch + FASTA(n).Multiplicity;
        Correctness.CorrectByMultiplicityExclFullMatch = Correctness.CorrectByMultiplicityExclFullMatch + FASTA(n).Multiplicity * score;
      end
      Correctness.CorrectSequences = Correctness.CorrectSequences + score;

      ROC(ed+1,1) = ROC(ed+1,1) + Counter * score;  % record degree of success
      ROC(ed+1,2) = ROC(ed+1,2) + Counter;          % count this sequence by its edit distance to its own group

      maxscore = max(GroupData(mn).OwnScore);

      if Verbose == 1,
        fprintf('Better: %3d Equal: %3d ',numbetter,numequal);
        fprintf('Score %4.2f', score);
        fprintf(' %20s', FASTA(n).Sequence);
        fprintf(' (%4d)', FASTA(n).Multiplicity);
        fprintf(' MLPS %6.2f', MLPS(n,mn,1));
        fprintf(' deficit %6.2f', maxscore-MLPS(n,mn,1));    % deficit in default order
        fprintf(' prct %6.2f', 100*Percentile(n,mn,1));      % percentile in default order
        fprintf(' CutScore %6.2f; ', CutoffScore(n,mn,1));
        fprintf(' Ed %2d,%2d', FullEditDistance(n,mn,1), CoreEditDistance(n,mn,1));  % edit distances in default order

        Text{seqnum} = [num2str(seqnum) char(9) num2str(2*gg) char(9)];
        Text{seqnum} = [Text{seqnum} GroupData(mn).MotifID char(9)];
        Text{seqnum} = [Text{seqnum} num2str(GroupData(mn).NumNT) char(9)];
        Text{seqnum} = [Text{seqnum} num2str(GroupData(mn).NumInstances) char(9)];
        Text{seqnum} = [Text{seqnum} num2str(GroupData(mn).NumBasepairs) char(9)];
        Text{seqnum} = [Text{seqnum} num2str(GroupData(mn).NumBPh) char(9)];
        Text{seqnum} = [Text{seqnum} num2str(GroupData(mn).NumBR) char(9)];
        Text{seqnum} = [Text{seqnum} num2str(GroupData(mn).Structured) char(9) char(9)];

        li = strfind(FASTA(n).Header,loopType);
        sp = strfind(FASTA(n).Header,' ');

        Text{seqnum} = [Text{seqnum} FASTA(n).Header(li(2):(li(2)+9)) char(9)];
        Text{seqnum} = [Text{seqnum} FASTA(n).Sequence char(9)];

        Text{seqnum} = [Text{seqnum} 'http://rna.bgsu.edu/rna3dhub/motif/view/' GroupData(mn).MotifID char(9)];
        Text{seqnum} = [Text{seqnum} GroupData(mn).Signature{1} char(9)];
        Text{seqnum} = [Text{seqnum} sprintf('%0.2f', MLPS(n,mn,1)) char(9)];
        Text{seqnum} = [Text{seqnum} sprintf('%0.2f', 100*Percentile(n,mn,r(mn))) char(9)];
        Text{seqnum} = [Text{seqnum} num2str(CoreEditDistance(n,mn,r(mn))) char(9)];

        Text{seqnum} = [Text{seqnum} num2str(numbetter) char(9)];

        if numbetter > 0,
          fprintf(' scores better against %3d groups: ', numbetter);
          for b = 1:min(NumBetterModelsToShow,numbetter),
            bmn = bettermodel(b);
            [thismlps,thisr] = max(Score(n,bmn,:));
            CoreEditDist = CoreEditDistance(n,bmn,thisr);
            FullEditDist = FullEditDistance(n,bmn,thisr);

            fprintf('%s, %2d NTs, ', GroupData(bmn).MotifID, GroupData(bmn).NumNT);
            fprintf('%-40s, Ed %2d,%2d, MLPS %5.2f, ', GroupData(bmn).Signature{1}, FullEditDist, CoreEditDist, MLPS(n,bmn,thisr));
            fprintf('deficit %5.2f, ', max(GroupData(bmn).OwnScore)-MLPS(n,bmn,thisr));
            fprintf('prct %6.2f; ', 100*Percentile(n,bmn,thisr));

            Text{seqnum} = [Text{seqnum} char(9) char(9) GroupData(bmn).MotifID char(9)];
            Text{seqnum} = [Text{seqnum} 'http://rna.bgsu.edu/rna3dhub/motif/view/' GroupData(bmn).MotifID char(9)];
            Text{seqnum} = [Text{seqnum} GroupData(bmn).Signature{1} char(9)];
            Text{seqnum} = [Text{seqnum} sprintf('%0.2f', MLPS(n,bmn,thisr)) char(9)];
            Text{seqnum} = [Text{seqnum} sprintf('%0.2f', 100*Percentile(n,bmn,thisr)) char(9)];
            Text{seqnum} = [Text{seqnum} num2str(CoreEditDist) char(9)];
            Text{seqnum} = [Text{seqnum} num2str(GroupData(bmn).NumNT) char(9)];
          end
          bestmodels(k) = bettermodel(1);
        elseif OwnCutoffMet,
          fprintf(' matches the original group, %s', GroupData(mn).Signature{1});
          bestmodels(k) = mn;
          Text{seqnum} = [Text{seqnum} 'OK'];
        else
          fprintf(' has no match');
          bestmodels(k) = 0;
        end

        fprintf('\n');

        % --------------- examine sequences with core edit distance 1 that are done right by edit distance but wrong by alignment score

        if 0 > 1 && Criterion == 2 && CoreEditDistance(n,mn,1) == 1 && OwnCutoffMet && numbetter == 0 && max(max(MLPS(n,:,:))) > MLPS(n,mn,1),
          fprintf('Core edit distance 1 success, but alignment score fail\n');
        end

        seqnum = seqnum + 1;
      end
    end

    if 0 > 1 && Verbose == 1 && length(bestmodels) > 0,
      [b,t,i] = zUniqueRows(bestmodels');
      fprintf('\n');
      for k = 1:length(b),
        if b(k) > 0,
          fprintf('Model %s was the best-scoring model for %4d sequences (%4.1f%%) %s', GroupData(str(b(k))).MotifID, t(k), 100*t(k)/sum(t), GroupData(str(b(k))).Signature{1});
          if b(k) == mn,
            fprintf(' <-- correct model');
          end
        else
          fprintf('There was no match for                          %4d sequences (%4.1f%%)', t(k), 100*t(k)/sum(t));
        end
        fprintf('\n');
      end
      fprintf('\n');
    end



  end
end

if Params.NumSequencesToShow < Inf,
  pause
end

Correctness.PercentageCorrectSequences = Correctness.CorrectSequences / Correctness.NumSeqs;
Correctness.PercentageCorrectMultiplicity = Correctness.CorrectByMultiplicity / Correctness.TotalMultiplicity;
Correctness.PercentageCorrectExclFullMatch = Correctness.CorrectByMultiplicityExclFullMatch / Correctness.TotalMultiplicityExclFullMatch;
Correctness.PercentageCorrectExclInteriorMatch = Correctness.CorrectByMultiplicityExclInteriorMatch / Correctness.TotalMultiplicityExclInteriorMatch;
Correctness.PercentageAcceptedMultiplicity = Correctness.AcceptedByMultiplicity / Correctness.TotalMultiplicity;
Correctness.PercentageAcceptedMultiplicityExclFullMatch = Correctness.AcceptedByMultiplicityExclFullMatch / Correctness.TotalMultiplicityExclFullMatch;
Correctness.PercentageAcceptedMultiplicityExclInteriorMatch = Correctness.AcceptedByMultiplicityExclInteriorMatch / Correctness.TotalMultiplicityExclInteriorMatch;
