% Makes a model based on Search and returns filenames of the
% model file, sequence file, and distribution file
% loopType must be specified ('IL' or 'HL')

function [Search,Node] = pMakeSingleJAR3DModel(Search,Param,Prior,loopType)

if nargin < 2,
  Param = [0 2 0 4 100 1 1 1 Inf 5 3];                   % See below

  % Parameters stored in Param:
  % Param(1) verbose
  % Param(2) method to use for basepair isostericity
  % Param(3) recognize extensible helices and model them as such
  % Param(4) adjust substitution probabilities for long-range interactions
  % Param(5) how far to look ahead for local basepair interactions
  % Param(6) use near interactions
  % Param(7) treat insertions as conserved bases
  % Param(8) normalize scores for insertions and basepairs:
  % Param(9) controls balance between isodiscrepancy and counts; L/(L+Param(9)) is the parameter for counts
  % Param(10) tells whether to use BPh interactions and the strength (recommend 5)
  % Param(11) tells whether to use BR interactions and the strenghth (recommend 3loa)
end

Verbose = Param(1);

if nargin < 3,
  Prior = [.5 .5 .5 .5 0];              % Prior distribution for insertion bases
end

% ----- Calculate coplanar measure for each File in Search

[L,N] = size(Search.Candidates);        % L = num instances; N = num NT
N = N - 1;                              % number of nucleotides

if ~isfield(Search.File(1),'Coplanar'),
    clear NewFile
    for ff = 1:length(Search.File),
        F = Search.File(ff);
        if ~isempty(F.NT),
            F.Coplanar = sparse(F.NumNT,F.NumNT);
            NewFile(ff) = F;
        end
    end
    Search.File = NewFile;

    for c = 1:length(Search.Candidates(:,1)),
        ff = Search.Candidates(c,N+1);
        i = Search.Candidates(c,1:N);

        for a = 1:length(i),
            for b = (a+1):length(i),
                if Search.File(ff).Edge(i(a),i(b)) ~= 0,
                    NT1 = Search.File(ff).NT(i(a));
                    NT2 = Search.File(ff).NT(i(b));
                    Pair = zAnalyzePair(NT1,NT2,CL,Exemplar);
                    Search.File(ff).Coplanar(i(a),i(b)) = Pair.Coplanar;
                    Search.File(ff).Coplanar(i(b),i(a)) = Pair.Coplanar;
                end
            end
        end
    end
end

% ----- Calculate loose coplanar measure for each File in Search

if ~isfield(Search.File(1),'LooseCoplanar'),
    clear NewFile
    for ff = 1:length(Search.File),
        F = Search.File(ff);
        if ~isempty(F.NT),
            F.LooseCoplanar = sparse(F.NumNT,F.NumNT);
            NewFile(ff) = F;
        end
    end
    Search.File = NewFile;

    for c = 1:length(Search.Candidates(:,1)),
        ff = Search.Candidates(c,N+1);
        i = Search.Candidates(c,1:N);

        for a = 1:length(i),
            for b = (a+1):length(i),
                if Search.File(ff).Edge(i(a),i(b)) ~= 0,
                    NT1 = Search.File(ff).NT(i(a));
                    NT2 = Search.File(ff).NT(i(b));
                    Pair = zLooseCoplanar(NT1,NT2,CL,Exemplar);
                    Search.File(ff).LooseCoplanar(i(a),i(b)) = Pair.Coplanar;
                    Search.File(ff).LooseCoplanar(i(b),i(a)) = Pair.Coplanar;
                end
            end
        end
    end
end

% --------------------------------------- Make SCFG/MRF model

[Node,Search] = pMakeMotifModelFromSSF(Search,Param,Prior,loopType,1:L);

if length(Node) > 0,
  OK = 1;
  for n = 1:length(Node),
    if strcmp(loopType,'IL') && strcmp(Node(n).type,'Junction'),
      OK = 0;
    end
  end
  if OK == 0,
    Node = [];
  end
end
