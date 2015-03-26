% pGenerateRandomMotifSequences generates random sequences to score against models

function [void] = pGenerateRandomMotifSequences(OutputBase,Release,Mode)

if nargin < 3,
  Mode = 3;                           % no constraints on base combinations
end

numrepeats = 50;                      % number of repetitions of each combination of strand lengths

loopType = Release(1:2);

switch loopType,
case 'HL'
  numfiles = 20;
case 'IL'
  numfiles = 20;
end

switch loopType,
case 'JL'
  Rotations = 3;                      % three rotations, for 3-way junctions
case 'IL'
  Rotations = 2;                      % two rotations are computed
case 'HL'
  Rotations = 1;                      % only one "rotation"
end

Release = strrep(Release,'\',filesep);
Release = strrep(Release,'/',filesep);
OutputPath = [OutputBase filesep Release];

ModelPath = [OutputPath filesep 'lib'];
SequencePath = ModelPath;

% ---------------------------------------- Read data from models

GroupData = pGetModelData(OutputPath,loopType);

clear SeqNames

for m = 1:length(GroupData),
  SeqNames{m} = [GroupData(m).MotifID '.fasta'];
end

[FASTA, OwnMotif] = pConcatenateFASTASequences(SequencePath, SeqNames);

TransitionFile = [OutputPath filesep 'transitions.mat'];

clear seqlengths
clear nummodelsbylengths

switch loopType,
case 'HL'
	for s = 1:length(FASTA),
		seqlengths(s,1) = length(FASTA(s).Sequence);
	end
  [uniqueseqlengths,i,j] = unique(seqlengths,'rows');
	for k = 1:length(i),
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
case 'IL'
	for s = 1:length(FASTA),
		i = strfind(FASTA(s).Sequence,'*');
		j = length(FASTA(s).Sequence) - i;
		seqlengths(s,1) = min(i-1,j);
		seqlengths(s,2) = max(i-1,j);
	end
  [uniqueseqlengths,i,j] = unique(seqlengths,'rows');

	for k = 1:length(i),
		a = uniqueseqlengths(k,1);
		b = uniqueseqlengths(k,2);
		nummodelsbylengths(a,b) = length(unique(OwnMotif(find(j==k))));
	end

  newlengths = [];

  for k = 1:length(uniqueseqlengths(:,1)),
    if uniqueseqlengths(k,1) < uniqueseqlengths(k,2),
      newlengths(end+1,:) = uniqueseqlengths(k,:) + [1 0];
    end
    newlengths(end+1,:) = uniqueseqlengths(k,:) + [0 1];
    newlengths(end+1,:) = uniqueseqlengths(k,:) + [1 1];
  end

  [uniqueseqlengths,i,j] = unique([uniqueseqlengths; newlengths],'rows');

	T = [seqlengths OwnMotif];
	T = sortrows(T,[1 2 3]);
end

% j maps from all sequences to rows of the uniqueseqlengths vector

fprintf('Unique sequence lengths:\n');
uniqueseqlengths

fprintf('Size of unique seq lengths:\n');
size(uniqueseqlengths)

seqlengths = [];

for r = 1:length(uniqueseqlengths(:,1)),
  seqlengths = [seqlengths; ones(numrepeats,1)*uniqueseqlengths(r,:)];
end

for v = 1:numfiles,
	Sequences = pMakeRandomSequencesWeighted(seqlengths,loopType,length(seqlengths(:,1)),TransitionFile,Mode);
  switch Mode
  case 1
    FN = [ModelPath filesep loopType '_RandomMotifSequencesNonCan_' num2str(v) '.fasta'];
  case 2
    FN = [ModelPath filesep loopType '_RandomMotifSequencesUAOK_' num2str(v) '.fasta'];
  case 3
    FN = [ModelPath filesep loopType '_RandomMotifSequences_' num2str(v) '.fasta'];
  end
	fid = fopen(FN,'w');
	for s = 1:length(Sequences),
		fprintf(fid,'%s\n',Sequences{s});
	end
	fclose(fid);
end


