% pUpdateModelWithSSF(Node,Search,...) uses the sequence variability in Search to adjust insertion and substitution probabilities

function [Node,Search] = pUpdateModelWithSSF(Node,Search,f,F,Param,Prior,loopType,File,UseCandidate,Normalize)

DeletionProbPerNucleotide = 0.005;     % deletion is less likely than putting a mismatched pair
InsertionProbOutsideFlank = 0.0001;

OutsideFlankLengthDist = [1-InsertionProbOutsideFlank InsertionProbOutsideFlank];

Verbose = Param(1);

[L,N] = size(Search.Candidates);        % L = num instances; N = num NT
N = N - 1;                              % number of nucleotides

Noncanonical = [];
switch loopType,
case 'HL'
  if N >= 4,
    Noncanonical = [2 N-1];
  end
case 'IL'
  if Search.Truncate > 2 && N - Search.Truncate > 2,  % number of nucleotides on each strand
    Noncanonical = [[2 N-1]; [Search.Truncate-1 Search.Truncate+2]];
  end
end

if nargin < 9,
    UseCandidate = 1:L;                     % Use all instances
else
    L = length(UseCandidate);               % Only use given instances
end

% ---------------------------- Set parameters for the nodes from instances

for n = 1:length(Node),
    switch Node(n).type
        case 'Initial' % ==========================================================
            if n == 1,
                % harsher penalties for first node insertions
                Node(n).leftLengthDist = OutsideFlankLengthDist;
                Node(n).rightLengthDist = OutsideFlankLengthDist;
            else
                a  = max(Node(n-1).LeftIndex);      % left NT of previous node
                aa = min(Node(n).LeftIndex);        % left NT of this node
                b  = min(Node(n-1).RightIndex);     % right NT of previous node
                bb = max(Node(n).RightIndex);       % right NT of this node

                % adjust insertion probs for initial nodes after the first node
                % ----------------------------- tally insertions on the left

                if n < length(Node)-1,              % no insertions after last pair
                  if a < aa,                        % possible insertion here

                    [LengthDist,LetterDist] = pInsertionDistribution(Search,a,UseCandidate,Prior,Normalize);
                    Node(n).leftLengthDist = LengthDist;
                    Node(n).leftLetterDist = LetterDist;

                  else
                    Node(n).leftLengthDist = [1];   % no insertions possible
                  end

                  if b > bb,
                    % ----------------------------- tally insertions on the right

                    [LengthDist,LetterDist] = pInsertionDistribution(Search,bb,UseCandidate,Prior,Normalize);
                    Node(n).rightLengthDist = LengthDist;
                    Node(n).rightLetterDist = LetterDist;

                  else
                    Node(n).rightLengthDist = [1];   % no insertions possible
                  end
                end
            end

        case 'Fixed' % ==========================================================

            Node(n).Delete = DeletionProbPerNucleotide;      % deletion probability controlled by number of nucleotides

            % ----------------------------- fixed base on the left
            letter = zeros(size(Prior));                     % record which bases occur

            if length(Node(n).leftLengthDist) > 1,           % fixed base on left
                a = Node(n).LeftIndex;                       % index within file
                for c = 1:L,                                 % loop through candidates
                    ff = Search.Candidates(UseCandidate(c),N+1); % file number
                    d = Search.Candidates(UseCandidate(c),a);    % index of interacting base
                    insb = Search.File(ff).NT(d).Code;       % A=1, C=2, G=3, U=4
                    if ~isempty(insb)
                        letter(insb) = letter(insb) + 1;
                    end
                end

                Pr = Prior;

                for e = 1:3,                                        % loop through edges
                    if Param(10) > 0 && sum(Search.BPh(a,:)==e) > 0,
                        Pr = Pr / Param(10);                  % weaken the prior distribution
%                        fprintf('pUpdateModelWithSSF: Fixed node gets new distribution on the left from BPh made by base %d using edge %d\n',a,e);
                    elseif Param(11) > 0 && sum(Search.BR(a,:)==e) > 0,
                        Pr = Pr / Param(11);                  % weaken the prior distribution
%                        fprintf('pUpdateModelWithSSF: Fixed node gets new distribution on the left from BR made by base %d using edge %d\n',a,e);
                    end
                end

                temp = letter + Pr;
                Node(n).leftLetterDist = temp / sum(temp);          % normalize

            end

            % ----------------------------- fixed base on the right
            letter = zeros(size(Prior));                     % record which bases occur

            if length(Node(n).rightLengthDist) > 1,    % fixed base on right
                a = Node(n).RightIndex;
                for c = 1:L,                               % loop through candidates
                    ff = Search.Candidates(UseCandidate(c),N+1); % file number
                    d = Search.Candidates(UseCandidate(c),a);    % index of interacting base
                    insb = Search.File(ff).NT(d).Code;       % A=1, C=2, G=3, U=4
                    if ~isempty(insb)
                        letter(insb) = letter(insb) + 1;
                    end
                end

                Pr = Prior;

                for e = 1:3,                                    % loop through edges
                    if Param(10) > 0 && sum(Search.BPh(a,:)==e) > 0,
                        Pr = Pr / Param(10);                  % weaken the prior distribution
                    elseif Param(11) > 0 && sum(Search.BR(a,:)==e) > 0,
                        Pr = Pr / Param(11);                  % weaken the prior distribution
                    end
                end

                temp = letter + Pr;
                Node(n).rightLetterDist = temp / sum(temp);          % normalize

            end

        case 'Basepair' % =======================================================

            AvgNumNT = 2;

            a = Node(n).LeftIndex;                   % left NT of the query motif
            b = Node(n).RightIndex;                  % right NT of the query motif

            if Verbose > 0,
                disp('pUpdateModelWithSSF: Getting consensus for a basepair');
            end

            Score = pConsensusPairSubstitution(a,b,f,Search.File,F,Search,Param,Normalize,Noncanonical);

            if Verbose > 0,
                fprintf('Original substitution probabilities\n');
                Node(n).SubsProb
                fprintf('Consensus substitution probabilities\n');
                Score
            end

            Node(n).SubsProb = Score;

            Search.SubsProb{a,b} = Score;
            Search.SubsProb{b,a} = Score';

            % ----------------------------- tally insertions on the left
            inscount = zeros(1,L);
            letter = Prior;                      % record which bases occur
            aa = min(Node(n+1).LeftIndex);       % next interacting base in the model
            bb = max(Node(n+1).RightIndex);      % next interacting base in the model
            if (n < length(Node)-1 || strcmp(loopType,'HL')) && aa <= bb,                    % no insertions after last pair

                % ----------------------------- tally insertions on the left
                [LengthDist,LetterDist] = pInsertionDistribution(Search,a,UseCandidate,Prior,Normalize);
                Node(n).leftLengthDist = LengthDist;
                Node(n).leftLetterDist = LetterDist;

                % ----------------------------- tally insertions on the right

if bb >= b,
  fprintf('pUpdateModelWithSSF:  Basepair and next node out of order! ************************ \n');
end
                [LengthDist,LetterDist] = pInsertionDistribution(Search,bb,UseCandidate,Prior,Normalize);
                Node(n).rightLengthDist = LengthDist;
                Node(n).rightLetterDist = LetterDist;

            end

            if n == length(Node)-1 && strcmp(loopType,'IL'), % last basepair in IL
                Node(n).leftLengthDist  = OutsideFlankLengthDist;
                Node(n).rightLengthDist = OutsideFlankLengthDist;
            end

            AvgNumNT = AvgNumNT + (0:(length(Node(n).leftLengthDist)-1)) * Node(n).leftLengthDist';
            AvgNumNT = AvgNumNT + (0:(length(Node(n).rightLengthDist)-1)) * Node(n).rightLengthDist';

            Node(n).Delete = DeletionProbPerNucleotide^AvgNumNT;              % deletion probability controlled by number of nucleotides

        case 'Cluster' % ===========================================================

            Indices = [Node(n).LeftIndex(Node(n).Left) Node(n).RightIndex(Node(n).Right)];
            FixedIndices = [Node(n).LeftIndex Node(n).RightIndex]; % indices of nucleotides not in basepairs
            Positions = 1:length(Indices);                         % positions within the cluster node
            AvgNumNT = length(Indices);

            % -------------- add 4x4 "interaction" matrices for fixed bases not in basepairs and not already accounted for; there should be none

            for ii = 1:length(Node(n).IBases(:,1)),
                a = Node(n).InterIndices(ii,1);
                b = Node(n).InterIndices(ii,2);

                FixedIndices = setdiff(FixedIndices,a);
                FixedIndices = setdiff(FixedIndices,b);
            end

            for fi = FixedIndices,
                [s,t] = size(Node(n).IBases);
                newinter = s + 1;                                % place to add new "interaction"
                dist = Prior;
                p = find(Indices == fi);                         % position within node
                Node(n).IBases(newinter,:) = [p p];              % record position within node
                Node(n).InterIndices(newinter,:) = [fi fi];      % index of this nucleotide in the motif
                Node(n).SubsProb(length(dist),length(dist),newinter) = 0; % make space if needed
                Node(n).SubsProb(:,:,newinter) = diag(dist);     % put dist down the diagonal
                Node(n).InteractionComment{newinter} = [' // Cluster node fixed base, nucleotide ' num2str(fi) ' at position ' num2str(p)];
                fprintf('pUpdateModelWithSSF:  Added a fixed base that was not previously modeled\n');
                pause
            end

            % ---------- find consensus 4x4 substitution matrix for basepairs and for fixed but non-basepairing nucleotides

            for ii = 1:length(Node(n).IBases(:,1)),
                a = Node(n).InterIndices(ii,1);
                b = Node(n).InterIndices(ii,2);

                if a ~= 0 && b ~= 0,
                    if Verbose > 0,
                        disp(['pUpdateModelWithSSF: Getting consensus for pairs in a cluster, bases ' num2str(a) ' and ' num2str(b)]);
                    end

                    Score = pConsensusPairSubstitution(a,b,f,Search.File,F,Search,Param,Normalize,Noncanonical,Prior);

                    if Verbose > 0,
                        fprintf('Original substitution probabilities\n');
                        Node(n).SubsProb(:,:,ii)

                        fprintf('Consensus substitution probabilities\n');
                        Score
                        fprintf('\n');
                    end

                    Node(n).SubsProb(:,:,ii) = Score;
                    Search.SubsProb{a,b} = Score;
                    Search.SubsProb{b,a} = Score';
                end
            end

            % ------- Adjust insertion probabilities between fixed nucleotides

            InternalNumbering = 1:(length(Node(n).LeftIndex) + length(Node(n).RightIndex));
            InsertionPositions = setdiff(InternalNumbering,length(Node(n).LeftIndex));
            InsertionPositions = InsertionPositions(1:(end-1));

            insertions = 0;                        % number of insertions identified so far

            for i = InsertionPositions,
                a = Indices(i);                    % index within the motif

                [LengthDist,LetterDist,insnumber] = pInsertionDistribution(Search,a,UseCandidate,Prior,Normalize);

                if max(insnumber) > 0,             % only put insertions where one is found
                    insertions = insertions + 1;
                    Node(n).Insertion(insertions).Position = i;
                    Node(n).Insertion(insertions).Index = a;
                    Node(n).Insertion(insertions).LengthDist = LengthDist;
                    Node(n).Insertion(insertions).LetterDist = LetterDist;
                    Node(n).InsertionComment{insertions} = ['// Insertion after nucleotide ' num2str(a) ', position ' num2str(i) ' in cluster node'];

                    AvgNumNT = AvgNumNT + (0:(length(LengthDist)-1)) * LengthDist';
                end
            end

            Node(n).Delete = DeletionProbPerNucleotide^AvgNumNT;              % deletion probability controlled by number of nucleotides

            Node(n).NormCons = pClusterNorm(Node(n).InterIndices,Node(n).SubsProb,Node(n).LeftIndex,Node(n).RightIndex);

        case 'Junction' % =========================================================

        case 'Hairpin'  % =========================================================

            if strcmp(loopType,'HL')         % only for "real" hairpins, not * hairpins

                a = Node(n).LeftIndex;
                b = Node(n).RightIndex;
                Positions = 1:(b-a+1);        % positions in the node, internal to the node
                Indices = a:b;                % indices of nucleotides in the motif
                FixedIndices = a:b;           % identify nts that are not in basepairs

                if b < a,                     % not supposed to happen
                    fprintf('pUpdateModelWithSSF:  Hairpin left and right indices out of order! **************\n');
                    a = min(a,b);             % should avoid crashing, but model won't be good
                    b = max(a,b);
                end

                % ---------- find consensus 4x4 substitution matrix for basepairs
                if ~isempty(Node(n).IBases),
                    [s,t] = size(Node(n).IBases);

                    for ii = 1:s,
                        if Node(n).IBases(ii,1) ~= 0 && Node(n).IBases(ii,2) ~= 0,
                            a = Node(n).InterIndices(ii,1);
                            b = Node(n).InterIndices(ii,2);

                            FixedIndices = setdiff(FixedIndices,a);
                            FixedIndices = setdiff(FixedIndices,b);

                            if Verbose > 0,
                                disp(['pUpdateModelWithSSF: Getting consensus for pairs in a cluster, bases ' num2str(a) ' and ' num2str(b)]);
                            end

                            Score = pConsensusPairSubstitution(a,b,f,Search.File,F,Search,Param,Normalize,Noncanonical);

                            if Verbose > 0,
                                fprintf('Original substitution probabilities\n');
                                Node(n).SubsProb(:,:,ii)

                                fprintf('Consensus substitution probabilities\n');
                                Score
                                fprintf('\n');
                            end

                            Node(n).SubsProb(:,:,ii) = Score;
                            Search.SubsProb{a,b} = Score;
                            Search.SubsProb{b,a} = Score';
                        end
                    end
                end

                % -------------- add 4x4 "interaction" matrices for fixed bases not in basepairs

                for fi = FixedIndices,
                    [s,t] = size(Node(n).IBases);
                    newinter = s + 1;                                % place to add new "interaction"
                    dist = Prior;
                    for j = 1:L                                      % loop through instances
                        f = Search.Candidates(UseCandidate(j),end);  % file number
                        k = Search.Candidates(UseCandidate(j),fi);   % nucleotide index in file
                        baseCode = Search.File(f).NT(k).Code;
                        dist(baseCode) = dist(baseCode) + 1;
                    end
                    dist = dist / sum(dist);                         % normalize
                    p = find(Indices == fi);                         % position within node
                    Node(n).IBases(newinter,:) = [p p];              % numbering internal to the node
                    Node(n).InterIndices(newinter,:) = [fi fi];      % index of this nucleotide in the motif
                    Node(n).SubsProb(length(dist),length(dist),newinter) = 0; % make space if needed
                    Node(n).SubsProb(:,:,newinter) = diag(dist);     % put dist down the diagonal
                    Node(n).InteractionComment{newinter} = [' // Hairpin fixed base, nucleotide ' num2str(fi) ' at position ' num2str(p)];
                end

                % -------------- check for insertions between fixed bases

                a = Node(n).LeftIndex;
                b = Node(n).RightIndex;
                insertions = 0;
                for i = a:(b-1),
                    [LengthDist,LetterDist,insnumber] = pInsertionDistribution(Search,i,UseCandidate,Prior,Normalize);

                    if max(insnumber) > 0,                                     % only insert if at least one insertion is found
                        insertions = insertions + 1;
                        Node(n).Insertion(insertions).Position = i-a+1;        % position within the hairpin
                        Node(n).Insertion(insertions).Index = i;               % index within the hairpin
                        Node(n).Insertion(insertions).LengthDist = LengthDist;
                        Node(n).Insertion(insertions).LetterDist = LetterDist;
                        Node(n).InsertionComment{insertions} = ['// Insertion in Hairpin after nucleotide ' num2str(i) ', position ' num2str(i-a+1)];
                    end
                end
            end
        end
    end
end
