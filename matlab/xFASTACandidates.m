% xFASTACandidates(File,Search,Direction,ModelName) shows a multiple sequence
% alignment of the candidates in Search.  Bases which correspond in the
% geometric search are aligned with one another.  If a maximum distance
% has been specified, or the maximum gap is small, bases between the
% bases in the candidate are also displayed.

% Direction can be +1 or -1; it tells the order in which to put the
% nucleotides in the first candidate.  0 means not to reorder them.

% Also returns alignment text in AText

function [Text,AText,BText,CText,FASTA,DText] = xFASTACandidates(File,Search,Direction,ModelName)

if nargin < 4,
  MN = Search.Query.Name;
else
  MN = ModelName;
end

AModelName = ModelName(1:6);       % temporary fix for filenames

Query      = Search.Query;
Candidates = Search.Candidates;

[L,N] = size(Candidates);
N = N - 1;

if Direction ~= 0,
  [y,p] = sort(Direction*double(Candidates(1,1:N)));
                                    % put nucleotides in inc/decreasing order
else
  p = 1:N;                          % do not re-order
end

Cand = double(Candidates(:,p));   % re-order nucleotides

F    = Candidates(:,N+1);           % file numbers

if isfield(Query,'MaxDiffMat'),
  MaxDiff = diag(Query.MaxDiffMat(p,p),1);
else
  MaxDiff = ones(1,N-1);
end

if Direction >= 0,
  MaxDiff = MaxDiff;
else
  MaxDiff = fliplr(MaxDiff);
end

% ---------------------------- Calculate maximum gaps between cand. nucleotides

maxinsert = zeros(1,N-1);
for c = 1:L,
  maxinsert = max(maxinsert,abs(diff(double(Cand(c,1:N))))-1);
end

% ---------------------------- Check about Truncate

if ~isfield(Search,'Truncate'),
  Search.Truncate = [];
else
  maxinsert = 0*maxinsert;
end

% ---------------------------- Print header line

t = 0;                                          % line of text we're on
r = 0;                                          % line of alignment we're on
k = N;                                          % number of columns being shown

% ----------------------------- Print alignment

CorrCodeList = zeros(c,N);                        % for bases corresp to query
for j = 1:N,
  InsCode{j}.Count = zeros(1,4);
  InsLength{j}.Count = [];
end

for c = 1:L,                                      % loop through candidates
  f = F(c);                                       % file number
  FN = strrep(Search.File(f).PDBFilename,'.pdb','');

  t = t + 1;
  Text{t} = [sprintf('> %s %5s %s%5s %s%5s', MN, File(f).Filename, File(f).NT(Cand(c,1)).Base, File(f).NT(Cand(c,1)).Number, File(f).NT(Cand(c,end)).Base, File(f).NT(Cand(c,end)).Number)];
  t = t + 1;
  Text{t} = '';                                   % start line of bases

  k = 1;                                          % column counter

  DText{c} = [sprintf('%s_Instance_%d has_name %s_%s_%s%s_%s%s', MN, c, MN, File(f).Filename, File(f).NT(Cand(c,1)).Base, File(f).NT(Cand(c,1)).Number, File(f).NT(Cand(c,end)).Base, File(f).NT(Cand(c,end)).Number)];

  for n = 1:N,                                    % print alignment for cand

    % ------------------------------------- add bases at conserved positions

    Text{t} = [Text{t} sprintf('%s', File(F(c)).NT(Cand(c,n)).Base)];
    j = File(F(c)).NT(Cand(c,n)).Code;
    CorrCodeList(c,n) = j;                        % store code of corresp base

    NT = Search.File(f).NT(Search.Candidates(c,n));
    MNum = NT.ModelNum;
    if isempty(MNum),
      MNum = 1;
    end
    r = r + 1;
    AText{r} = sprintf('%s_Instance_%d_Column_%d_%s corresponds_to_sequence %s_Sequence_%d_Position_%d_%s', ModelName, c, n, NT.Base, ModelName, c, length(Text{t}), NT.Base);
    % Note:  The nucleotide IDs in the following line might not be formed correctly, especially for nucleotides from symmtery operations
    BText{r} = sprintf('%s_Instance_%d_Column_%d_%s corresponds_to_PDB %s|%d|%s|%s|%s', ModelName, c, n, NT.Base, FN, MNum, NT.Chain, NT.Number, NT.Base);
    CText{r} = sprintf('%s_Instance_%d_Column_%d_%s corresponds_to_group %s_Column_%d', ModelName, c, n, NT.Base, ModelName, n);

    % ----------------------- add inserted bases between conserved positions

    if n < N,
      k = k + 1;
      if (MaxDiff(n) == Inf) || (maxinsert(n) >= 4) || any(n+1 == Search.Truncate),
        Text{t} = [Text{t} sprintf('*')];            % show break in strand
      else
        if Cand(c,n+1) - Cand(c,n) > 1,             % increasing order
          ic = 1;                                   % insertion counter
          for i = (Cand(c,n)+1):(Cand(c,n+1)-1),
            Text{t} = [Text{t} sprintf('%c', File(F(c)).NT(i).Base)];   % show insertions
            j = File(F(c)).NT(i).Code;
            InsCode{n}.Count(j) = InsCode{n}.Count(j) + 1;

            NT = Search.File(f).NT(i);
            MNum = NT.ModelNum;
            if isempty(MNum),
              MNum = 1;
            end
            r = r + 1;
            AText{r} = sprintf('%s_Instance_%d_Column_%d-%d_Insertion_%d_%s corresponds_to_sequence %s_Sequence_%d_Position_%d_%s', ModelName, c, n, n+1, ic, NT.Base, ModelName, c, length(Text{t}), NT.Base);
            if isfield(NT,'ID') && ~isempty(NT.ID),
              BText{r} = sprintf('%s_Instance_%d_Column_%d-%d_Insertion_%d_%s corresponds_to_PDB %s_%d_%s_%s_%s', ModelName, c, n, n+1, ic, NT.Base, FN, MNum, NT.Chain, NT.Number, NT.Base);
            else
              BText{r} = sprintf('%s_Instance_%d_Column_%d-%d_Insertion_%d_%s corresponds_to_PDB %s_%d_%s_%s_%s', ModelName, c, n, n+1, ic, NT.Base, FN, MNum, NT.Chain, NT.Number, NT.Base);
            end
            CText{r} = sprintf('%s_Instance_%d_Column_%d-%d_Insertion_%d_%s corresponds_to_group %s_Column_%d-%d_Insertion', ModelName, c, n, n+1, ic, NT.Base, ModelName, n, n+1);
            ic = ic + 1;                            % increment insertion count
          end
        elseif Cand(c,n+1) - Cand(c,n) < -1,        % decreasing order
          ic = 1;                                   % insertion counter
          for i = (Cand(c,n)-1):-1:(Cand(c,n+1)+1),
            Text{t} = [Text{t} sprintf('%c', File(F(c)).NT(i).Base)];   % show insertions
            j = File(F(c)).NT(i).Code;
            InsCode{n}.Count(j) = InsCode{n}.Count(j) + 1;

            NT = Search.File(f).NT(i);
            MNum = NT.ModelNum;
            if isempty(MNum),
              MNum = 1;
            end
            r = r + 1;
            AText{r} = sprintf('%s_Instance_%d_Column_%d-%d_Insertion_%d_%s corresponds_to_sequence %s_Sequence_%d_Position_%d_%s', ModelName, c, n, n+1, ic, NT.Base, ModelName, c, length(Text{t}), NT.Base);
            if isfield(NT,'ID') && ~isempty(NT.ID),
              BText{r} = sprintf('%s_Instance_%d_Column_%d-%d_Insertion_%d_%s corresponds_to_PDB %s_%d_%s_%s_%s', ModelName, c, n, n+1, ic, NT.Base, FN, MNum, NT.Chain, NT.Number, NT.Base);
            else
              BText{r} = sprintf('%s_Instance_%d_Column_%d-%d_Insertion_%d_%s corresponds_to_PDB %s_%d_%s_%s_%s', ModelName, c, n, n+1, ic, NT.Base, FN, MNum, NT.Chain, NT.Number, NT.Base);
            end
            CText{r} = sprintf('%s_Instance_%d_Column_%d-%d_Insertion_%d_%s corresponds_to_group %s_Column_%d-%d_Insertion', ModelName, c, n, n+1, ic, NT.Base, ModelName, n, n+1);
            ic = ic + 1;                            % increment insertion count
          end
        end

        for i=1:(1 + maxinsert(n) - abs(Cand(c,n+1)-Cand(c,n))),
          Text{t} = [Text{t} sprintf('-')];
          k = k + 1;
        end

      end

      h = abs(Cand(c,n+1)-Cand(c,n));                % number of insertions + 1
      if h > length(InsLength{n}.Count),
        InsLength{n}.Count(h) = 1;
      else
        InsLength{n}.Count(h) = InsLength{n}.Count(h) + 1;
      end
    end
  end
%  Text{t} = [Text{t} sprintf('%s', File(F(c)).NT(Cand(c,N)).Base)];
%  j = File(F(c)).NT(Cand(c,N)).Code;
%  CorrCodeList(c,N) = j;                           % store last corresp base


  c = 1;
  n = 1;
  while c < length(Text),
    FASTA(n).Header = Text{c};
    c = c + 1;
    FASTA(n).Sequence = Text{c};
    c = c + 1;
    n = n + 1;
  end

  drawnow

end
