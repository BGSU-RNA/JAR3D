% pAddPairstoFASTAforDiagnostic adds the specified pairs of letters to a collection of sequences in the data structure FASTA

function [FASTA] = pAddPairstoFASTAforDiagnostic(FASTA,NumSequences,DiagnosticMode)

if any(DiagnosticMode == [2 4]),
  for n = 1:NumSequences,
    s = FASTA(n).Sequence;
    i = strfind(s,'*');             % locations of * characters
    for j = length(i):-1:1,
      s = [s(1:(i(j)-2)) 'U' s((i(j)-1):(i(j)+1)) 'U' s((i(j)+2):end)];
    end

fprintf('pJAR3DDiagnostics %20s becomes %20s\n', FASTA(n).Sequence, s);

    FASTA(n).Sequence = s;
  end
end

if any(DiagnosticMode == [3 4]),
  for n = 1:NumSequences,
    s = FASTA(n).Sequence;
    s = [s(1) 'U' s(2:(end-1)) 'U' s(end)];

fprintf('pJAR3DDiagnostics %20s becomes %20s\n', FASTA(n).Sequence, s);

    FASTA(n).Sequence = s;
  end
end
