% pGenerateRandomMotifSequences generates random sequences to score against models

function [numfiles] = pGenerateRandomMotifSequences(OutputBase,Release,Mode)

if nargin < 3
  Mode = 3;                           % no constraints on base combinations
end

loopType = Release(1:2);

switch loopType
case 'HL'
  numrepeats = 50;                      % number of repetitions of each combination of strand lengths
  numfiles = 20;
case 'IL'
  numrepeats = 50;                      % number of repetitions of each combination of strand lengths
  numfiles = 20;
case 'J3'
  numrepeats = 30;
  numfiles = 20;
case 'J4'
  numrepeats = 25;
  numfiles = 20;
otherwise
  numrepeats = 25;
  numfiles = 20;
end

Release = strrep(Release,'\',filesep);
Release = strrep(Release,'/',filesep);
OutputPath = [OutputBase filesep Release];

ModelPath = [OutputPath filesep 'lib'];
SequencePath = ModelPath;

% ---------------------------------------- Read data from models

GroupData = pGetModelData(OutputPath,loopType);

clear SeqNames

for m = 1:length(GroupData)
  SeqNames{m} = [GroupData(m).MotifID '.fasta'];
end

[FASTA, OwnMotif] = pConcatenateFASTASequences(SequencePath, SeqNames);

TransitionFile = [OutputPath filesep 'transitions.mat'];

clear seqlengths
clear nummodelsbylengths

switch loopType
case 'HL'
	for s = 1:length(FASTA)
		seqlengths(s,1) = length(FASTA(s).Sequence);
	end
  [uniqueseqlengths,i,j] = unique(seqlengths,'rows');
	for k = 1:length(i)
		a = uniqueseqlengths(k,1);
		nummodelsbylengths(a,1) = length(unique(OwnMotif(find(j==k))));
	end

  newlengths = [];

  for k = 1:length(uniqueseqlengths(:,1)),
    newlengths(end+1,:) = uniqueseqlengths(k,:) + 1;
  end

  [uniqueseqlengths,i,j] = unique([uniqueseqlengths; newlengths],'rows');

	T = [seqlengths OwnMotif];
	T = sortrows(T,[1 2]);

otherwise
  % IL, J3, J4, etc.
  for s = 1:length(FASTA)
		i = strfind(['*' FASTA(s).Sequence '*'],'*');
    lengths = diff(i) - 1;
    % put shortest length first by permuting cyclically
    [y,i] = min(lengths);
    i = i(1);
    seqlengths(s,:) = lengths([i:end 1:(i-1)]);
    % seqlengths(s,1) = min(i-1,j);
		% seqlengths(s,2) = max(i-1,j);
	end
  [uniqueseqlengths,i,j] = unique(seqlengths,'rows');

  % not used
	% for k = 1:length(i)
	% 	a = uniqueseqlengths(k,1);
	% 	b = uniqueseqlengths(k,2);
	% 	nummodelsbylengths(a,b) = length(unique(OwnMotif(find(j==k))));
	% end

  % round out the selection of sequences somewhat
  newlengths = [];
  if strcmp(loopType,'IL')
    for k = 1:length(uniqueseqlengths(:,1))
      if uniqueseqlengths(k,1) < uniqueseqlengths(k,2)
        newlengths(end+1,:) = uniqueseqlengths(k,:) + [1 0];
      end
      newlengths(end+1,:) = uniqueseqlengths(k,:) + [0 1];
      newlengths(end+1,:) = uniqueseqlengths(k,:) + [1 1];
    end
  else
    for k = 1:length(uniqueseqlengths(:,1))
      for m = 1:length(uniqueseqlengths(k,:))
        t = uniqueseqlengths(k,:);
        t(m) = t(m) + 1;
        newlengths(end+1,:) = t;
      end
    end
  end

  [uniqueseqlengths,i,j] = unique([uniqueseqlengths; newlengths],'rows');

	T = [seqlengths OwnMotif];
	T = sortrows(T,[1 2 3]);
end

% j maps from all sequences to rows of the uniqueseqlengths vector

fprintf('pGenerateRandomMotifSequences: Unique sequence lengths:\n');
disp(uniqueseqlengths)

fprintf('pGenerateRandomMotifSequences: Size of unique seq lengths:\n');
size(uniqueseqlengths)

seqlengths = [];

for r = 1:length(uniqueseqlengths(:,1))
  seqlengths = [seqlengths; ones(numrepeats,1)*uniqueseqlengths(r,:)];
end

for v = 1:numfiles
  switch Mode
  case 1
    FN = [ModelPath filesep loopType '_RandomMotifSequencesNonCan_' num2str(v) '.fasta'];
  case 2
    FN = [ModelPath filesep loopType '_RandomMotifSequencesUAOK_' num2str(v) '.fasta'];
  case 3
    FN = [ModelPath filesep loopType '_RandomMotifSequences_' num2str(v) '.fasta'];
  end

  % only make and write the file if it does not already exist;
  % avoid overwriting previously generated ones
  if ~isfile(FN)
    Sequences = pMakeRandomSequencesWeighted(seqlengths,loopType,length(seqlengths(:,1)),TransitionFile,Mode);
    fid = fopen(FN,'w');
    for s = 1:length(Sequences)
      fprintf(fid,'%s\n',Sequences{s});
    end
    fclose(fid);
  end
end
