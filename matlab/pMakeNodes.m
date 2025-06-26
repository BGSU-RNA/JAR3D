% pMakeNodes(File,NTNumber,LastNTNumber,Truncate,Interact,Node,n) makes a secondary
% structure node model based on the Edge interaction matrix in File, starting at
% NTNumber and ending at LastNTNumber.  It assigns various nodes consistent with
% this secondary structure.
% Truncate indicates the first nucleotide of each new strand, for modeling
% internal and junction loops.
% Interact, Node, and n are optional parameters specified when pMakeNodes is called by itself,
% for example, after a Junction is found.

function [Node] = pMakeNodes(File,Param,NTNumber,LastNTNumber,Truncate,Data,Node,n)

if nargin < 2
    Verbose = 1;
end

if nargin < 3
    NTNumber = 1;
end

Verbose         = Param(1);
% method          = 4;             % method for assigning pair subst probs
Extension       = 1;             % whether to extend stems with no LR inter
% AdjustSubsForLR = 1;             % adjust ins, basepair subs probs for LR inter
cdepth          = 20;            % how far to look ahead for a cluster; can be large for motif atlas
usenear         = 0;             % use near pairs for basepair probabilities
% insertionconserved = 0;          % treat insertions as conserved bases?
% jcdepth         = 4;             % how far to look for a junction cluster
% Normalize = 1;

% Parameters stored in Param:
% Param(1) verbose
% Param(2) method to use for basepair isostericity
% Param(3) recognize extensible helices and model them as such
% Param(4) adjust substitution probabilities for long-range interactions
% Param(5) how far to look ahead for local basepair interactions
% Param(6) use near interactions
% Param(7) treat insertions as conserved bases
% Param(8) normalize scores for insertions and basepairs

if length(Param) > 6
    Normalize = Param(8);
end

if length(Param) > 6
    insertionconserved = Param(7);
end

if length(Param) > 5
    usenear = Param(6);
end

if length(Param) > 4
    cdepth = Param(5);
end

if length(Param) > 3
    AdjustSubsForLR = Param(4);
end

if length(Param) > 2
    Extension = Param(3);
end

if length(Param) > 1
    method  = Param(2);
end

if nargin < 5
    Truncate = [];
end

if nargin < 8
    n=0;                            % current node number
end

% temporary to show lots of information
Verbose = 1;


% -------------------------- if File is a text string (filename), load the file

if strcmp(class(File),'char')
    Filename = File;
    File = zGetNTData(Filename,0);
end

% ----------------- if NTNumber is a cell array of numbers, look up the indices

if strcmp(class(NTNumber),'char')
    NTNumber = {NTNumber};
end

if strcmp(class(NTNumber),'cell')
    NTNumber = zIndexLookup(File,NTNumber);
end

% ------------------------------------------ Set key variables

N = length(File.NT);                       % number of nucleotides in File
DelProb = 0.01;                            % nominal deletion probability, used in other programs
% for basepairs
TertiaryFreeNode = 0;                      % first node in this stem making
% no tertiary interactions beyond it

if ~isfield(File,'BasePhosphate')
    File.BasePhosphate = sparse(zeros(N,N));
end

if ~isfield(File,'BaseRibose')
    File.BaseRibose = sparse(zeros(N,N));
end

if nargin < 4
    LastNTNumber = N;
elseif strcmp(class(LastNTNumber),'cell')
    LastNTNumber = zIndexLookup(File,LastNTNumber);
elseif strcmp(class(LastNTNumber),'char')
    LastNTNumber = zIndexLookup(File,{LastNTNumber});
end

% ------------------------------------------ Store indices of interacting bases

load PairExemplars                              % used in other programs

if nargin < 6

    E = abs(fix(File.Edge));                    % don't distinguish subcategories
    G = E .* (E < 16) .* (E ~= 0);              % consider basepairing only,
                                                % including bifurcated, water-ins

    CanonicalPairs = {'CG','GC','AU','UA','GU','UG'};
    CanonicalcWW = sparse(zeros(N,N));

    [i,j,k] = find(G == 1);                     % cWW basepairs
    for ii = 1:length(i),
        if ismember([File.NT(i(ii)).Base File.NT(j(ii)).Base],CanonicalPairs),
            CanonicalcWW(i(ii),j(ii)) = 1;
            CanonicalcWW(j(ii),i(ii)) = 1;
        end
    end

    if usenear > 0
        G = G +  E .* (E > 100) .* (E < 116);   % near basepairs too
    end

    H = (G ~= 0) .* max(File.Crossing == 0, abs(G) == 1);
    % 1 for nested pairs, 0 otherwise

    J = abs(G .* (File.Crossing >  0));         % long-range basepairs only

    GG = G;
    for i = 1:(length(GG(:,1))-1)
        if pStrandBreaksBetween(i,i+1,Truncate) == 0
            GG(i,i+1) = 0;                      % eliminate pairs btw adjacent
            GG(i+1,i) = 0;
        end
    end

    for a = 1:N                                 % loop through nucleotides
        k = find(G(a,:));                       % find indices of interacting bases
        [y,L] = sort(E(a,k));                   % sort by edge interaction category
        Interact{a}.Categ = abs(File.Edge(a,k(L)));   % store categories
        Interact{a}.Index = k(L);               % store indices of interacting bases
    end

    % ------------------------------------------ prepare to identify motifs

    HasMotif     = zeros(1,length(File.NT));
    HasGUPacking = zeros(1,length(File.NT));   % highly-conserved motif

    if isfield(File,'Nucl')
        for i = 1:length(File.NT)
            if ~isempty(File.Nucl(i).Motif)
                HasMotif(i) = 1;
                for m = 1:length(File.Nucl(i).Motif)
                    if ~isempty(strfind(File.Nucl(i).Motif(m).Name,'GU_packing'))
                        HasGUPacking(i) = 1;
                    end
                end
            end
        end
    end

    Data.Interact = Interact;
    Data.HasMotif = HasMotif;
    Data.HasGUPacking = HasGUPacking;
    Data.E = E;
    Data.G = G;
    Data.H = H;
    Data.J = J;
    Data.CanonicalcWW = CanonicalcWW;
else
    Interact = Data.Interact;
    HasMotif = Data.HasMotif;
    HasGUPacking = Data.HasGUPacking;
    E            = Data.E;
    G            = Data.G;
    H            = Data.H;
    J            = Data.J;
    CanonicalcWW = Data.CanonicalcWW;
end

% ------------------------------------------ Set up initial values of counters
a  = NTNumber;                             % first index; current index
A  = a;                                    % previous interacting base on left
AA = a;                                    % previous cWW base on left

B  = LastNTNumber;                         % next base on right
BB = a;                                    % previous cWW base on right

if Verbose > 0
    fprintf('pMakeNodes starting loop:     a=%s B=%s\n', File.NT(a).ID, File.NT(B).ID);
end

% Initial node creation -------------------------------------------------

n = n+1;                                   % move to next node

pMakeNodesNewNode;                         % set up blank node with all fields
Node(n).type      = 'Initial';             % node type
Node(n).nextnode  = n+1;                   % index of next node in tree
Node(n).LeftIndex = a;                     % index of first base on left
Node(n).RightIndex= B;                     % index of last base on right
Node(n).Comment = ' // Initial node';

if G(a,B) == 0
    pMakeNodesProbeForInsertions;              % probe for insertions, each strand
                                               % moves a and B to next basepair
end

% ---------------------------------------------------------------------------

EndLoop = 0;                               % flag for the end of the loop

while (EndLoop == 0) && (a <= LastNTNumber) % while not the end of the loop,

    % Strategy for junctions: push all the way to where there are no more basepairs to be modeled
    % on the current stem, because maybe there will be ones to model on the next stem

    StartJunction = 0;
    if pStrandBreaksBetween(a,B,Truncate) > 1
        % there is a junction ahead, is now the time to start it?

        % indices of last nucleotides on current left and right strands
        aaa = pNextTruncation(a,Truncate,N);
        BBB = pNextTruncation(B,Truncate,1);

        % nucleotides on next strand(s) (NSL and NSR are the same for J3)
        [NSL,NSR] = pNextStrands(a,B,Truncate,N);

        if a == aaa
            % reached the end of this strand
            StartJunction = 1;
        elseif sum(sum(H(a,NSL))) > sum(sum(H(a,BBB:B)))
            % more interactions with NSL than right strand, stop
            StartJunction = 1;
        end

        if B == BBB
            % reached the end of this strand
            StartJunction = 1;
        elseif sum(sum(H(B,NSR))) > sum(sum(H(B,a:aaa)))
            % more interactions with NSL than right strand, stop
            StartJunction = 1;
        end
    end

    if StartJunction                         % determined by pMakeNodesProbeForInsertions

        pMakeNodesJunction                   % make model for junction

        % fprintf('Inserted a junction, pausing before moving on to the next motif\n')
        % pause

        return                               % nothing left to do, junctions cover the rest

    else                                     % not a junction

        % ---------------------------------- Identify basepair or cluster

        [aaa,BBB] = pNextNucleotides(a,B,Truncate,N,cdepth);

        LS = (a+1):aaa;                         % left strand after a
        RS = BBB:(B-1);                         % right strand before B

        if Verbose > 10
            fprintf('pMakeNodes left strand can go up    as far as %s\n', File.NT(aaa).ID);
            fprintf('pMakeNodes right strand can go back as far as %s\n', File.NT(BBB).ID);
        end

        if HasMotif(a)   % ---------------------- Insert motif model when needed

            pMakeNodesMotif

        elseif (H(a,B) > 0) && pStrandBreaksBetween(a,B,Truncate) == 1 && ...
            (a == pNextTruncation(a,Truncate,N)) && (B == pNextTruncation(B,Truncate,1))
            % a and B pair, a is at the end of a strand, and B is at the start of the next strand

            % fprintf('pMakeNodes making a basepair:  a=%s B=%s\n', File.NT(a).ID, File.NT(B).ID);
            pMakeNodesBasepair                       % add basepair with insertions
            % fprintf('pMakeNodes after  a basepair:  a=%s B=%s\n', File.NT(a).ID, File.NT(B).ID);

        elseif (H(a,B) > 0) && ((B - a) > 1 || pStrandBreaksBetween(a,B,Truncate) > 0) && ...
            ((cdepth == 0) || (sum(sum(G(a,[LS RS]))) == 0 && sum(sum(G([LS RS],B))) == 0))
            % a and B interact, but not also with other nearby bases

            % fprintf('pMakeNodes making   basepair:  a=%s B=%s\n', File.NT(a).ID, File.NT(B).ID);
            pMakeNodesBasepair                       % add basepair with insertions
            % fprintf('pMakeNodes after    basepair:  a=%s B=%s\n', File.NT(a).ID, File.NT(B).ID);

        else     % a and B also interact with nearby bases - use a cluster node

            % [a aaa BBB B]
            % zShowInteractionTable(File,unique([a:aaa BBB:B]));

            if pStrandBreaksBetween(a,B,Truncate) > 0
                fprintf('pMakeNodes making a cluster near strand break a=%s B=%s\n', File.NT(a).ID, File.NT(B).ID);
                pMakeNodesCluster;       % across a strand break; internal loop maybe with complicated interactions
            else
                fprintf('pMakeNodes checking on making a cluster a=%s B=%s\n', File.NT(a).ID, File.NT(B).ID);
                pMakeNodesCluster;

                if make_hairpin == 1
                    % the cluster would cover all nucleotides in the current stem, so use a hairpin instead
                    if Verbose > 0
                        fprintf('pMakeNodes making a hairpin a=%s B=%s\n', File.NT(a).ID, File.NT(B).ID);
                    end
                    pMakeNodesHairpin;
                    EndLoop = 1;
                end
            end
        end                                          % basepair or cluster

        % ------------------- check for tertiary interactions in this stem

        pMakeNodesCheckForExtensibility
        if Verbose > 10
            fprintf('pMakeNodes checked extensibi:  a=%s B=%s\n', File.NT(a).ID, File.NT(B).ID);
        end

        % ------------------- check for truncation and hairpin
        GG = G;
        for i = 1:(length(GG(:,1))-1)
            if pStrandBreaksBetween(i,i+1,Truncate) == 0
                GG(i,i+1) = 0;                   % eliminate pairs between adjacent, like cSH
                GG(i+1,i) = 0;
            end
        end

        % -------- Have not added a hairpin but maybe we need to stop anyway
        if (EndLoop == 0)
            % removed ismember(a-1,Truncate) || but added a > B for after a cluster
            if (a >= B && length(Truncate) > 0) || ismember(a,Truncate) || isempty(File.NT(a).Base)

                % fprintf('pMakeNodes truncating with * HL:  a=%s B=%s\n', File.NT(a).ID, File.NT(B).ID);

                pMakeNodesTruncate                % add * hairpin

            elseif (a >= B) || ((sum(sum(abs(GG(a:B,a:B)))) == 0))
                % no nucleotides left, or no further basepairs except cSH
                if (TertiaryFreeNode > 0) && isempty(Truncate) && Extension > 0
                    % extensible region, not a truncated model, only applies to HL
                    if Verbose > 0
                        fprintf('pMakeNodes extra basepairs: a=%s B=%s\n', File.NT(a).ID, File.NT(B).ID);
                    end
                        pMakeNodesExtraBasepairs;       % add extra basepairs if extensible
                end

                if Verbose > 0
                    fprintf('pMakeNodes adding hairpin:     a=%s B=%s\n', File.NT(a).ID, File.NT(B).ID);
                end

                pMakeNodesHairpin

            elseif G(a,B) > 0
                fprintf('No need to probe because of interaction between a=%s B=%s\n', File.NT(a).ID, File.NT(B).ID);
            else

                fprintf('pMakeNodes probe for inserts:  a=%s B=%s\n', File.NT(a).ID, File.NT(B).ID);
                pMakeNodesProbeForInsertions   % add insertions to basepair node or new initial node if needed
                fprintf('pMakeNodes probed for insert:  a=%s B=%s\n', File.NT(a).ID, File.NT(B).ID);

            end                                  % hairpin or insertions
        elseif TertiaryFreeNode > 0 && isempty(Truncate) && Extension > 0
            pMakeNodesExtraBasepairs;     % add extra basepairs if extensible, only for HL
        end                               % if EndLoop == 0
    end                                   % if junction or not
end                                       % while (EndLoop == 0) & (a <= N),

% ---------------------------------- Poisson distribution for lengths -------

function [d] = subPoisson(m)

n = max(3,2*m);

d = exp(-m) * (m .^ (0:n)) ./ factorial(0:n);

d = d / sum(d);                     % normalize
