% pGetModelData loads data on JAR3D models

function [GroupData, MotifEquivalence] = pGetModelData(OutputPath,loopType)

fprintf('pGetModelData\n');
filename = [OutputPath filesep loopType '_GroupData.mat'];
load(filename,'GroupData');
fprintf('Loaded GroupData from %s\n',filename);

try
  filename = [OutputPath filesep loopType '_GroupData_with_full_cutoffs.mat'];
  load(filename,'GroupData');
  fprintf('Loaded GroupData with model-specific cutoffs from %s\n',filename);
end

keep = ones(1,length(GroupData));
for m = 1:length(GroupData),
  if isempty(GroupData(m).MotifID),
    keep(m) = 0;
  end
end
GroupData = GroupData(find(keep));

MotifEquivalence = zeros(length(GroupData),length(GroupData));
MotifNotEquiv    = zeros(length(GroupData),length(GroupData));
FN = [OutputPath filesep 'lib' filesep 'equivalent_motifs.txt'];

switch loopType,
case 'IL'
  Rotations = 1;
case 'HL'
  Rotations = 0;
end

for m = 1:length(GroupData),
  MotifNames{m,1} = GroupData(m).MotifID;
  SequenceLengths = [];
  for s = 1:length(GroupData(m).OwnSequence),
    SequenceLengths(s) = length(GroupData(m).OwnSequence{s}) - Rotations;
  end
  GroupData(m).SequenceLengths = SequenceLengths;
  GroupData(m).MeanSequenceLength = mean(SequenceLengths);

  if strcmp(GroupData(m).Signature{1},'cWW-cWW') && strcmp(loopType,'IL'),
    GroupData(m).Structured = 0;
  elseif strcmp(GroupData(m).Signature{1},'cWW') && strcmp(loopType,'HL'),
    GroupData(m).Structured = 0;
  else
    GroupData(m).Structured = 1;
  end
%  fprintf('%d %s %s\n',GroupData(m).Structured,GroupData(m).MotifID,GroupData(m).Signature{1});

end

if exist(FN,'file') > 0,
  fid = fopen(FN,'r');
  L = '';
  while ischar(L),
    L = fgetl(fid);
    a = zStringSplit(L,char(9));
    if length(a) >= 3,
      i = find(ismember(MotifNames,a(1)));
      j = find(ismember(MotifNames,a(2)));
      if length(i) == 1 && length(j) == 1,
        if str2num(a{3}) > 0,
          MotifEquivalence(i,j) = 1;
          MotifEquivalence(j,i) = 1;
        else
          MotifNotEquiv(i,j) = 1;
          MotifNotEquiv(j,i) = 1;
        end
      elseif length(i) == 0,
        fprintf('Motif %s appears in the list of equivalences but not in the list of motifs\n',a{1});
      elseif length(j) == 0,
        fprintf('Motif %s appears in the list of equivalences but not in the list of motifs\n',a{2});
      end
    end
  end
  fclose(fid);
end

for i = 1:length(GroupData),
  MotifEquivalence(i,i) = 1;
end

MotifEquivalence = MotifEquivalence^20;
MotifEquivalence = (MotifEquivalence > 0);
MotifEquivalence = MotifEquivalence - MotifNotEquiv;  % put -1 where the file asserts non-equivalence

return

[NumModels,MotifNames,ModNames,ModNums] = pGetModelNames(ModelPath,loopType);

NumModels = length(MotifNames);

for m = 1:NumModels,
  GroupData(m).MotifID = MotifNames{m};

  fid = fopen([ModelPath filesep MotifNames{m} '_data.txt']);
  s = fgetl(fid);                 % read signature(s)
  i = strfind(s,' ');             % find spaces in signature, if any
  i = [1 i length(s)];
  for j = 1:(length(i)-1),
    GroupData(m).Signature{j} = s(i(j):i(j+1));
  end

  b = fgetl(fid);                 % read number of nucleotides
  b = regexprep(b,'[a-z,A-Z]','');
  GroupNumNT(m) = str2num(b);     % convert to integer
  GroupData(m).NumNT = str2num(b);
  a = fgetl(fid);                 % phonetic string

  a = fgetl(fid);
  a = regexprep(a,'[a-z,A-Z]','');
  GroupData(m).NumBasepairs = str2num(a);
  a = fgetl(fid);
  a = regexprep(a,'[a-z,A-Z]','');
  GroupData(m).NumStacks = str2num(a);
  a = fgetl(fid);
  a = regexprep(a,' base-phosphate','');
  GroupData(m).NumBPh = str2num(a);
  a = fgetl(fid);
  a = regexprep(a,' base-ribose','');
  GroupData(m).NumBR = str2num(a);
  a = fgetl(fid);
  a = regexprep(a,'[a-z,A-Z]','');
  GroupData(m).NumInstances = str2num(a);
  a = fgetl(fid);
  a = regexprep(a,'[a-z,A-Z]','');
  GroupData(m).Truncate = str2num(a);

  fclose(fid);

  fid = fopen([ModelPath filesep MotifNames{m} '_ownscores.txt']);
  s = fgetl(fid);
  c = 1;
  while ischar(s),
    [token, remain] = strtok(s,char(9));
    GroupData(m).OwnSequence{c} = token;
    GroupData(m).OwnScore(c) = str2num(remain);
    c = c + 1;
    s = fgetl(fid);
  end
  fclose(fid);

  if  GroupData(m).NumBPh > 0 || GroupData(m).NumBR > 0,
    Structured(m) = 1;
    GroupData(m).Structured = 1;
  elseif GroupData(m).NumBasepairs > 2 && strcmp(loopType,'IL'),
    Structured(m) = 1;
    GroupData(m).Structured = 1;
  elseif GroupData(m).NumBasepairs > 1 && strcmp(loopType,'HL'),
    Structured(m) = 1;
    GroupData(m).Structured = 1;
  else
    Structured(m) = 0;
    GroupData(m).Structured = 0;
  end

  si = strrep(GroupData(m).Signature{1},'L','R');
  GroupData(m).NumFixed = length(strfind(si,'-R-'));

end


if 0 > 1,
  [i,j] = find(MotifEquivalence);
  for k = 1:length(i),
    fprintf('%s = %s\n',MotifNames{i(k)},MotifNames{j(k)});
  end
end
