% Makes a model based on Search and returns filenames of the
% model file, sequence file, and distribution file
% loopType must be specified ('IL' or 'HL')

function [filePaths,Node] = pWriteSingleJAR3DModel(Search,origmatfilename,Verbose,NumRandomSeqs)

if nargin < 3,
  Verbose = 0;
end

if nargin < 4,
  NumRandomSeqs = 100000;
end

loopType = origmatfilename(1:2);

MN = strrep(origmatfilename,'.mat','');

Param = [1 2 0 1 100 1 1];                   % See below

% Parameters stored in Param:
% Param(1) verbose
% Param(2) method to use for basepair isostericity
% Param(3) recognize extensible helices and model them as such
% Param(4) adjust substitution probabilities for long-range interactions
% Param(5) how far to look ahead for local basepair interactions
% Param(6) use near interactions
% Param(7) treat insertions as conserved bases

Prior = [.5 .5 .5 .5 0];              % Prior distribution for insertion bases

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
    save(FN,'Search','-mat');
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
    save(FN,'Search','-mat');
end

% --------------------------------------- Notify

if Verbose > 0,
  fprintf('pMakeModelsFromLibrary:  Making a JAR3D SCFG/MRF model for %s\n', MN);
end

% --------------------------------------- Make model and write it

% try

[Node,Truncate,Signature,RSignature] = pMakeMotifModelFromSSF(Search,Param,Prior,loopType,1:L);

% catch
%   full(min(30,abs(Search.File.Edge)))
%   Signature = 'trouble';
% end

Search.Signature = Signature;
Search.RSignature = RSignature;           %
% save(FN,'Search','-mat');


if ~strcmp(MN(4:5),'99'),                 % don't do this for helices!
    for n = 1:length(Node),
        if fix(abs(Node(n).Edge)) == 1,       % cWW basepair
            %        Node(n).SubsProb = ones(4,4)/16;    % make cWW pairs noninformative
            %        Node(n).SubsProb = [0 0 0 1; 0 0 1 0; 0 1 0 1; 1 0 1 0]/6;
            % make cWW pairs noninformative
        end
    end
end

Search.Truncate = Truncate;

% --------------------------------------- Write sequences in FASTA format
Text = xFASTACandidates(Search.File,Search,0,MN);

fprintf('pMakeModelsFromLibrary:  First sequence is %s\n',Text{2});

fid = fopen(['Sequences' filesep MN '.fasta'],'w');
for t = 1:length(Text),
    fprintf(fid,'%s\n',Text{t});
end
fclose(fid);

Search.ownsequencefasta = [MN '.fasta'];

% --------------------------------------- Write model in text file

disp(['pWriteSingleJAR3DModel: MN is ' MN]);

Text = pNodeToSCFGModelText(Node,5);

% pJavaNodeFile(Search.Query,Node,5,['Models' filesep MN],0);

fid = fopen(['Models' filesep MN '_model.txt'],'w');
for i = 1:length(Text),
  fprintf(fid,'%s\n', Text{i});
end
fclose(fid);

Search.modelfilename = [MN '.txt'];
% save(FN,'Search','-mat');

% --------------------------------------- Generate and write empirical
% --------------------------------------- score distribution
top10 = 0;
sampSize = NumRandomSeqs;
precision = 4 ;           % ------------- Number of decimal places emp dist is precise to

%if sampSize > 0,
  D = pModelLengthDist(Node);
  if ~(exist(['Models' filesep 'Emp. Distributions']) == 7),
    mkdir(['Models' filesep 'Emp. Distributions']);
  end
  distFN = ['Models' filesep 'Emp. Distributions' filesep MN '.txt'];
  BetterEmpDist(MN,loopType,sampSize,precision,top10,distFN,D);
%end

filePaths.seq = [pwd filesep 'Sequences' filesep MN '.fasta'];
filePaths.mod = [pwd filesep 'Models' filesep MN '_model.txt'];
filePaths.dist = [pwd filesep distFN];

