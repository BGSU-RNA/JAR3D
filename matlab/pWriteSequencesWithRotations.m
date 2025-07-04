% pWriteSequencesWithRotations writes the sequences in FASTA to a file, and if they
% are from an internal loop, rotates the sequences and writes them again, etc.

function [AllSequencesFile,rotatedFASTA] = pWriteSequencesWithRotations(DiagnosticPath,FASTA,loopType,Rotations,CF)

% in the future, use loopType to determine the number of rotations,
% especially if you want to do every permutation of the strands for J3, J4, etc.

if nargin < 5
  CF = '';
end

r = 1;                          % write first rotation
rotatedFASTA{r} = FASTA;
AllSequencesFile{1} = [DiagnosticPath filesep CF '_All_Sequences_' num2str(r) '.fasta'];
fid = fopen(AllSequencesFile{1},'w');
for n = 1:length(FASTA)
  fprintf(fid,'>%s\n',FASTA(n).Header);
  fprintf(fid,'%s\n',FASTA(n).Sequence);
end
fclose(fid);

if Rotations > 1                       % more than one rotation, for IL, J3
  FASTA_R = FASTA;                     % start with the same FASTA data each time

  for r = 2:Rotations
    for n = 1:length(FASTA_R)
      a = FASTA_R(n).Sequence;
      i = strfind(a,'*');
      i = i(r-1);                      % pick one of the * characters to rotate around
      a = [a((i+1):end) '*' a(1:(i-1))];

      b = FASTA_R(n).Aligned;
      i = strfind(b,'*');
      i = i(r-1);                      % pick one of the * characters to rotate around
      b = [b((i+1):end) '*' b(1:(i-1))];

      FASTA_R(n).Sequence = a;
      FASTA_R(n).Aligned = b;
    end

    rotatedFASTA{r} = FASTA_R;

    AllSequencesFile{r} = [DiagnosticPath filesep CF '_All_Sequences_' num2str(r) '.fasta'];
    fid = fopen(AllSequencesFile{r},'w');
    for n = 1:length(FASTA_R)
      fprintf(fid,'>%s\n',FASTA_R(n).Header);
      fprintf(fid,'%s\n',FASTA_R(n).Sequence);
    end
    fclose(fid);
  end
end