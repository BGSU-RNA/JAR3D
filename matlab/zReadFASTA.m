% zReadFASTA reads the indicated fasta file and stores the records in an array.  Each element has three fields:
% Data(n).Header is the header for the nth sequence
% Data(n).Aligned is the sequence, with gaps
% Data(n).Sequence is the sequence with all gaps stripped out

function [Data] = zReadFASTA(Filename, RemoveJAR3DSymbols, Verbose, KeepSequence)

if nargin < 2,
  RemoveJAR3DSymbols = 0;
end

if nargin < 3,
  Verbose = 0;
end

if nargin < 4,
  KeepSequence = 1;
end

fid = fopen(Filename,'r');

if fid > 0

  L = 'a';
  c = 0;

  while ischar(L),
    L = fgetl(fid);
    if ischar(L),
      if isempty(L),
        L(1) = '-';
      end
      if L(1) == '>',
        c = c + 1;
        Data(c).Header = L(2:end);
        Data(c).Aligned = '';
      else
        Data(c).Aligned = [Data(c).Aligned L];
      end
    end
  end

  fclose(fid);

  for n = 1:length(Data),
    le(n) = length(Data(n).Aligned);
    Data(n).Aligned = strrep(Data(n).Aligned,'.','-'); % dots to hyphens
    Data(n).Aligned = strrep(Data(n).Aligned,'T','U'); % DNA to RNA
    Data(n).Aligned = strrep(Data(n).Aligned,' ','');  % remove spaces
    if RemoveJAR3DSymbols > 0, 
      Data(n).Aligned = strrep(Data(n).Aligned,'(','');  % remove other chars
      Data(n).Aligned = strrep(Data(n).Aligned,')','');  % remove other chars
      Data(n).Aligned = strrep(Data(n).Aligned,'{','');  % remove other chars
      Data(n).Aligned = strrep(Data(n).Aligned,'}','');  % remove other chars
      Data(n).Aligned = strrep(Data(n).Aligned,'<','');  % remove other chars
      Data(n).Aligned = strrep(Data(n).Aligned,'>','');  % remove other chars
      Data(n).Aligned = strrep(Data(n).Aligned,'|','');  % remove other chars
    end
%    Data(n).Aligned  = [Data(n).Aligned '-'];  % extra column for NT in struct
                                               % with nothing in FASTA
    if KeepSequence == 1,
      Data(n).Sequence = strrep(Data(n).Aligned,'-','');
    end
  end

  if min(le) < max(le),
    if Verbose > 0,
      disp('Sequences have different lengths')
      [min(le) max(le)]
    end
    for n = 1:length(Data),
      for c = (le(n)+1):max(le),
        Data(n).Aligned(c) = '-';
      end
    end
  end

  if Verbose > 0,
    for n = 1:length(Data),
      fprintf('>%s\n', Data(n).Header);
      fprintf('%s\n', Data(n).Sequence);
    end
  end  

else

  fprintf('zReadFASTA: Could not open file %s\n', Filename);

end
