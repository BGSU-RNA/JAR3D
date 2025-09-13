% pWriteSequencesWithRotations writes the sequences in FASTA to a file
% If they are from an IL, a separate file with the strands rotated is produced
% If from J3, a sequence like AA*BB*CC is rotation 0, BB*CC*AA is rotation 1, CC*AA*BB is rotation 2
% Similarly with J4, J5, etc.

function [AllSequencesFile,rotatedFASTA] = pWriteSequencesWithRotations(DiagnosticPath,FASTA,loopType,max_rotation,CF)

% in the future, use loopType to determine the number of rotations,
% especially if you want to do every permutation of the strands for J3, J4, etc.

if nargin < 5
  CF = '';
end

r = 0;                                 % write file with 0 rotation
rotatedFASTA{1} = FASTA;             % matlab indices start at 1
AllSequencesFile{1} = [DiagnosticPath filesep CF '_All_Sequences_' num2str(r) '.fasta'];
fid = fopen(AllSequencesFile{1},'w');
for n = 1:length(FASTA)
  fprintf(fid,'>%s\n',FASTA(n).Header);
  fprintf(fid,'%s\n',FASTA(n).Sequence);
end
fclose(fid);

n = 1;
a = FASTA(n).Sequence;
fprintf('max_rotation=%d r=%d n=%d a=%s FASTA(n).Sequence=%s\n',max_rotation,r,n,a,FASTA(n).Sequence);

if max_rotation > 0                    % more than one rotation, for IL, J3
  FASTA_R = FASTA;                     % start with the same FASTA data each time

  for r = 1:max_rotation
    for n = 1:length(FASTA_R)
      a = FASTA(n).Sequence;
      i = strfind(a,'*');
      i = i(r);                        % r tells which * in the sequence to cycle around
      a = [a((i+1):end) '*' a(1:(i-1))];

      if n == 1
        fprintf('max_rotation=%d r=%d n=%d a=%s FASTA(n).Sequence=%s\n',max_rotation,r,n,a,FASTA(n).Sequence);
      end

      b = FASTA(n).Aligned;
      i = strfind(b,'*');
      i = i(r);                        % r tells which * in the sequence to cycle around
      b = [b((i+1):end) '*' b(1:(i-1))];

      FASTA_R(n).Sequence = a;
      FASTA_R(n).Aligned = b;
    end

    rotatedFASTA{r+1} = FASTA_R;

    AllSequencesFile{r+1} = [DiagnosticPath filesep CF '_All_Sequences_' num2str(r) '.fasta'];
    fid = fopen(AllSequencesFile{r+1},'w');
    for n = 1:length(FASTA_R)
      fprintf(fid,'>%s\n',FASTA_R(n).Header);
      fprintf(fid,'%s\n',FASTA_R(n).Sequence);
    end
    fclose(fid);
  end
end