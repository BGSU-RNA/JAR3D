% pMakeRandomSequencesWeighted(seqlengths,loopType,sampsize,TransitionFile) generates random sequences

% BetterEmpDist('IL_225_05_cWW-cWW-cSH','IL',100)

function [Text,D,FASTA] = pMakeRandomSequencesWeighted(seqlengths,loopType,sampsize,TransitionFile,Mode)

if nargin < 5
  Mode = 3;                       % no constraints on base combinations
end

[s,t] = size(seqlengths);

if s == sampsize
    fprintf('pMakeRandomSequencesWeighted: using passed in sequence dimensions\n')
    dim = seqlengths';                  % first argument is the desired sequence lengths
else
    Node = seqlengths;            % no idea when that would be needed
    D = pModelLengthDist(Node);
    cWW_M = zeros(6,1);           %CG,GC,AU,UA,GU,UG
    initial_M = zeros(4,1);       %A,C,G,U
    transition_M = zeros(4,4);    %(A,C,G,U)^2

    % ---------------------------------------- generate strand lengths
    switch loopType
        case 'HL'
            % dim=zeros(1,sampsize);
            dim = myrandsample(D,1,sampsize)-1;
        case 'IL'
            ld = sum(D,2);                       % marginal distribution of left
            dim=zeros(2,sampsize);

            dim(1,:) = myrandsample(ld,sampsize,1) - 1;

            for i = 1:sampsize
                rd = D(dim(1,i)+1,:);            % almost the conditional distn
                rd = rd / sum(rd);               % normalize
                dim(2,i) = myrandsample(rd)-1;
            end
    end
end

fprintf('pMakeRandomSequencesWeighted: dim variable columns 1 to 10\n')
disp(dim(:,1:10))

load(TransitionFile)
fprintf('pMakeRandomSequencesWeighted: cWW distribution, initial letter distribution, transition matrix')
disp(cWW_M)
disp(initial_M)
disp(transition_M)

% -------------------------------------- make sequences

r = 1;
v = 1;
cWW = ['CG','GC','AU','UA','GU','UG'];          % string 'CGGCAUUAGUUG'
[strands,N]=size(dim);                          % N is the sample size
Text{N} = '';                                   % allocate space for all sequences
FASTA(N).Header = '';

switch Mode
case 1
    % prohibit these base combinations next to the flanking pair
    RemoveCombinations = ['CG'; 'GC'; 'AU'; 'UA'; 'GU'; 'UG'];
case 2
    RemoveCombinations = ['CG'; 'GC'; 'AU'; 'GU'; 'UG'];          % UA tWH occurs often
case 3
    RemoveCombinations = '';
end

letter = 'ACGU';

retry = 0;

if strands == 1
    % generate flanking cWW and initial nucleotide for all samples
    openprand = myrandsample(cWW_M,sampsize,1);
    mrand     = myrandsample(initial_M,sampsize,1);

    for i = 1:N
        Text{r} = sprintf('> Variant %d',i);
        r = r + 1;
        openp = openprand(i,1); % flanking cWW pair
        LL = max(dim(1,i),2);

        s = GenerateStrand(LL,mrand(i,1),transition_M);

        if LL > 2
            while ismember([s(1) s(end)],RemoveCombinations,'rows')
                retry = retry + 1;
                s = GenerateStrand(LL,mrand(i,1),transition_M);
            end
        end

        %fprintf('pMakeRandomSequencesWeighted: LL=%2d s=%s\n',LL,s);

        Text{r} = [cWW(2*openp-1) s cWW(2*openp)];
        r = r + 1;

        FASTA(v).Header = sprintf('> Variant %d',i);
        FASTA(v).Sequence = [cWW(2*openp-1) s cWW(2*openp)];
        v = v + 1;
    end
elseif strands == 2
    % IL
    openprand  = myrandsample(cWW_M,sampsize,1);
    closeprand = myrandsample(cWW_M,sampsize,1);

    for i = 1:N                                 % N sequences
        Text{r} = sprintf('> Variant %d',i);
        r = r + 1;

        LL = max(dim(1,i),2);                   % left strand length
        RL = max(dim(2,i),2);                   % right strand length
        openp = openprand(i,1);
        closep = closeprand(i,1);

        if LL == 2
            s = [];
        elseif LL == 3
            m = myrandsample(initial_M,1,1);
            s = letter(m);
        else
            m = myrandsample(initial_M,1,1);
            s = letter(m);
            for j = 1:(LL-3)
                m = myrandsample(transition_M(m,:));
                s = [s letter(m)];
            end
        end

        if RL == 2
            t = [];
        elseif RL == 3
            m = myrandsample(initial_M,1,1);
            t = letter(m);
        else
            m = myrandsample(initial_M,1,1);
            t = letter(m);
            for j = 1:(RL-3)
                m = myrandsample(transition_M(m,:));
                t = [t letter(m)];
            end
        end

        if LL > 2 && RL > 2
            while ismember([s(1) t(end)],RemoveCombinations,'rows') || ismember([t(1) s(end)],RemoveCombinations,'rows'),
                retry = retry + 1;
                if LL == 2
                    s = [];
                elseif LL == 3
                    m = myrandsample(initial_M,1,1);
                    s = letter(m);
                else
                    m = myrandsample(initial_M,1,1);
                    s = letter(m);
                    for j = 1:(LL-3),
                        m = myrandsample(transition_M(m,:));
                        s = [s letter(m)];
                    end
                end

                if RL == 2
                    t = [];
                elseif RL == 3
                    m = myrandsample(initial_M,1,1);
                    t = letter(m);
                else
                    m = myrandsample(initial_M,1,1);
                    t = letter(m);
                    for j = 1:(RL-3),
                        m = myrandsample(transition_M(m,:));
                        t = [t letter(m)];
                    end
                end
            end
        end

        Text{r} = [cWW(2*openp-1) s cWW(2*closep) '*' cWW(2*closep-1) t cWW(2*openp)];

% fprintf('%5d %2d %2d %s * %s | %s\n',i,LL,RL,s,t,Text{r});

        FASTA(v).Header = sprintf('> Variant %d',i);
        FASTA(v).Sequence = Text{r};
        v = v + 1;
        r = r + 1;

    end

else
    % J3, J4, and will also cover IL and maybe HL
    retry = 0;
    for i = 1:N                              % N sequences
        % generate and check interior sequences
        badPairFound = 1;                    % pretend, so we get into the while loop

        while badPairFound
            interior = cell(1,strands);
            for strand = 1:strands
                sL = max(dim(strand,i),2);        % strand length

                % generate the interior of each strand sequence
                if sL == 2
                    s = '';
                elseif sL == 3
                    m = myrandsample(initial_M,1,1);
                    s = letter(m);
                else
                    s = repmat('X',1,sL-2);
                    m = myrandsample(initial_M,1,1);
                    s(1) = letter(m);
                    for j = 2:(sL-2)
                        m = myrandsample(transition_M(m,:));
                        s(j) = letter(m);
                    end
                end
                interior{strand} = s;
            end
            interior{strands+1} = interior{1};  % simplify code below

            badPairFound = 0;

            % check to see if any of the pairs are "bad"
            % if so, set badPairFound = 1
            if strands == 1
                s = interior{1};
                if length(s) >= 2
                    if ismember([s(1) s(end)],RemoveCombinations,'rows')
                        badPairFound = 1;
                    end
                end
            else
                for strand = 1:strands
                    s = interior{strand};
                    t = interior{strand+1};  % wraps around to first strand
                    if length(s) > 0 && length(t) > 0
                        if ismember([s(end) t(1)],RemoveCombinations,'rows')
                            badPairFound = 1;
                        end
                    end
                end
            end

            retry = retry + badPairFound;
        end

        % now that we have a good sequence, record it
        f = '';                         % full sequence
        for strand = 2:strands
            % generate code for one flanking pair
            p = myrandsample(cWW_M,1,1);
            % add this pair; the 2*p-1 goes at the beginning of a strand
            f = [f cWW(2*p) '*' cWW(2*p-1) interior{strand} ]; %#ok<AGROW>
        end
        f = [cWW(2*p-1) interior{1} f cWW(2*p)]; %#ok<AGROW>

        Text{r} = sprintf('> Variant %d',i);
        r = r + 1;
        Text{r} = f;
        r = r + 1;

        FASTA(v).Header = sprintf('> Variant %d',i);
        FASTA(v).Sequence = f;
        v = v + 1;
    end
end

if Mode ~= 3
    fprintf('pMakeRandomSequencesWeighted: Generating %d sequences required %d additional attempts to avoid the specified pairs\n',N,retry);
end

function [s] = GenerateStrand(LL,m,transition_M)

    letter = 'ACGU';
    if LL == 2
        s = '';
    elseif LL == 3
        s = letter(m);
    else
        s = letter(m);
        for w = 1:(LL-3)
            m = myrandsample(transition_M(m,:));
            s = [s letter(m)];
        end
    end

