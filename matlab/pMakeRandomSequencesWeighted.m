% pMakeRandomSequencesWeighted(Node,loopType,sampsize,TransitionFile) generates random sequences

% BetterEmpDist('IL_225_05_cWW-cWW-cSH','IL',100)

function [Text,D,FASTA] = pMakeRandomSequencesWeighted(Node,loopType,sampsize,TransitionFile,AvoidCanonical)

if nargin < 5,
    AvoidCanonical = 0;
end

[s,t] = size(Node);

if s == sampsize,
    dim = Node';                   % first argument is the desired sequence lengths
else
    D = pModelLengthDist(Node);
    cWW_M = zeros(6,1);           %CG,GC,AU,UA,GU,UG
    initial_M = zeros(4,1);       %A,C,G,U
    transition_M = zeros(4,4);    %(A,C,G,U)^2

    % ---------------------------------------- generate strand lengths
    switch loopType,
    case 'IL'
        ld = sum(D,2);                       % marginal distribution of left
        dim=zeros(2,sampsize);

        dim(1,:) = myrandsample(ld,sampsize,1) - 1;

        for i = 1:sampsize,
            rd = D(dim(1,i)+1,:);            % almost the conditional distn
            rd = rd / sum(rd);               % normalize
            dim(2,i) = myrandsample(rd)-1;
        end
    case 'HL'
        dim=zeros(1,sampsize);
        dim = myrandsample(D,1,sampsize)-1;
    end
end

load(TransitionFile)

% -------------------------------------- make sequences

r = 1;
v = 1;
cWW = ['CG','GC','AU','UA','GU','UG'];
[strands,N]=size(dim);
Text{N} = '';                                   % allocate space for all sequences
FASTA(N).Header = '';

if AvoidCanonical > 0,
    Canonical = ['CG'; 'GC'; 'AU'; 'UA'; 'GU'; 'UG'];
else
    Canonical = '';
end

letter = 'ACGU';

retry = 0;

if strands == 2,
    openprand  = myrandsample(cWW_M,sampsize,1);
    closeprand = myrandsample(cWW_M,sampsize,1);

    for i = 1:N,                                % N sequences
        LL = max(dim(1,i),2);                   % left strand length
        RL = max(dim(2,i),2);                   % right strand length
        Text{r} = sprintf('> Variant %d',i);
        r = r + 1;
        openp = openprand(i,1);
        closep = closeprand(i,1);

        if LL == 2,
            s = [];
        elseif LL == 3,
            m = myrandsample(initial_M,1,1);
            s = letter(m);
        else,
            m = myrandsample(initial_M,1,1);
            s = letter(m);
            for j = 1:(LL-3),
                m = myrandsample(transition_M(m,:));
                s = [s letter(m)];
            end
        end

        if RL == 2,
            t = [];
        elseif RL == 3,
            m = myrandsample(initial_M,1,1);
            t = letter(m);
        else,
            m = myrandsample(initial_M,1,1);
            t = letter(m);
            for j = 1:(RL-3),
                m = myrandsample(transition_M(m,:));
                t = [t letter(m)];
            end
        end

        if LL > 2 && RL > 2,
            while ismember([s(1) t(end)],Canonical,'rows') || ismember([s(end) t(1)],Canonical,'rows'),
                retry = retry + 1;
                if LL == 2,
                    s = [];
                elseif LL == 3,
                    m = myrandsample(initial_M,1,1);
                    s = letter(m);
                else,
                    m = myrandsample(initial_M,1,1);
                    s = letter(m);
                    for j = 1:(LL-3),
                        m = myrandsample(transition_M(m,:));
                        s = [s letter(m)];
                    end
                end

                if RL == 2,
                    t = [];
                elseif RL == 3,
                    m = myrandsample(initial_M,1,1);
                    t = letter(m);
                else,
                    m = myrandsample(initial_M,1,1);
                    t = letter(m);
                    for j = 1:(RL-3),
                        m = myrandsample(transition_M(m,:));
                        t = [t letter(m)];
                    end
                end
            end
        end

        Text{r} = [cWW(2*openp-1) s cWW(2*closep-1) '*' cWW(2*closep) t cWW(2*openp)];

% fprintf('%5d %2d %2d %s * %s | %s\n',i,LL,RL,s,t,Text{r});

        FASTA(v).Header = sprintf('> Variant %d',i);
        FASTA(v).Sequence = Text{r};
        v = v + 1;

        r = r + 1;

    end
else
    openprand = myrandsample(cWW_M,sampsize,1);
    mrand     = myrandsample(initial_M,sampsize,1);

    for i = 1:N,
        Text{r} = sprintf('> Variant %d',i);
        r = r + 1;
        openp = openprand(i,1); % flanking cWW pair
        LL = max(dim(1,i),2);

        s = GenerateStrand(LL,mrand(i,1),transition_M);

        if LL > 2,
            while ismember([s(1) s(end)],Canonical,'rows'),
                retry = retry + 1;
                s = GenerateStrand(LL,mrand(i,1),transition_M);
            end
        end

%fprintf('%2d %s\n',LL,s);

        Text{r} = [cWW(2*openp-1) s cWW(2*openp)];
        r = r + 1;

        FASTA(v).Header = sprintf('> Variant %d',i);
        FASTA(v).Sequence = [cWW(2*openp-1) s cWW(2*openp)];
        v = v + 1;
    end
end

if AvoidCanonical > 0,
    fprintf('Generating %d sequences required %d additional attempts to avoid canonical pairs\n',N,retry);
end

function [s] = GenerateStrand(LL,m,transition_M);

    letter = 'ACGU';
    if LL == 2,
        s = [];
    elseif LL == 3,
        s = letter(m);
    else,
        s = letter(m);
        for i = 1:(LL-3),
            m = myrandsample(transition_M(m,:));
            s = [s letter(m)];
        end
    end
