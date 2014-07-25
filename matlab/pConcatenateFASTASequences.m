% pConcatenateFASTASequences reads files listed in SeqNames from SequencePath and stores them in a data structure called FASTA.  OwnGroup has the same dimensions as FASTA, and each entry tells the number of the entry in SeqNames to which it corresponds.

function [FASTA, OwnGroup] = pConcatenateFASTASequences(SequencePath, SeqNames)

for n = 1:length(SeqNames),
  newFASTA = zReadFASTA([SequencePath filesep SeqNames{n}]);
  m = length(newFASTA);
  g = n*ones(m,1);                % group that these sequences belong to
  if n == 1,                      % first set of sequences
    FASTA = newFASTA;
    OwnGroup = g;
  else
    FASTA = [FASTA newFASTA];     % concatenate fasta lines
    OwnGroup = [OwnGroup; g];     % mark each line by group it's from
  end
end
