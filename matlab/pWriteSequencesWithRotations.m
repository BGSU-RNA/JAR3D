% pWriteSequencesWithRotations writes the sequences in FASTA to a file, and if they are from an internal loop, rotates the sequences and writes them again, etc.

function [AllSequencesFile] = pWriteSequencesWithRotations(ModelPath,FASTA,loopType,Rotations,CF)

if nargin < 5,
  CF = '';
end

r = 1;                          % first rotation
AllSequencesFile{1} = [ModelPath filesep CF '_All_Sequences_' num2str(r) '.fasta'];
fid = fopen(AllSequencesFile{1},'w');
for n = 1:length(FASTA),
  fprintf(fid,'>%s\n',FASTA(n).Header);
  fprintf(fid,'%s\n',FASTA(n).Sequence);
end  
fclose(fid);

if Rotations > 1,                       % more than one rotation, for IL, JL
  FASTA_R = FASTA;
end

for r = 2:Rotations,
  for n = 1:length(FASTA_R),
    a = FASTA_R(n).Sequence;
    b = FASTA_R(n).Aligned;

    % fprintf('%s and %s become ', a, b);

    i = strfind(a,'*');
    a = [a((i(1)+1):end) '*' a(1:(i(1)-1))];
    i = strfind(b,'*');
    b = [b((i(1)+1):end) '*' b(1:(i(1)-1))];

    % fprintf('%s and %s.\n', a, b);

    FASTA_R(n).Sequence = a;
    FASTA_R(n).Aligned = b;
  end

  AllSequencesFile{r} = [ModelPath filesep CF '_All_Sequences_' num2str(r) '.fasta'];
  fid = fopen(AllSequencesFile{r},'w');
  for n = 1:length(FASTA_R),
    fprintf(fid,'>%s\n',FASTA_R(n).Header);
    fprintf(fid,'%s\n',FASTA_R(n).Sequence);
  end  
  fclose(fid);
end
