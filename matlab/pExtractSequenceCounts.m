function [cWWCounts, initialCounts, transitionCounts, coreCounts, initialTransitionCounts] = pExtractSequenceCounts(Text)
cWWCounts = zeros(6,1);           %CG,GC,AU,UA,GU,UG
initialCounts = zeros(4,1);       %A,C,G,U
coreCounts = zeros(4,1);          %A,C,G,U
transitionCounts = zeros(4,4);    %(A,C,G,U)^2
initialTransitionCounts = zeros(4,4);    %(A,C,G,U)^2
numSeq = floor(length(Text)/2);
for i = 1:numSeq,
    sequence = Text{2*i};
    seqLength = length(sequence);
    clear breakPoint;
    breakPoint = strfind(sequence, '*');
    cWW1 = strcat(sequence(1),sequence(seqLength));
    switch cWW1,
        case {'CG','cg','Cg','cG'}
            cWWCounts(1) = cWWCounts(1) + 1;
        case {'GC','gc','Gc','gC'}
            cWWCounts(2) = cWWCounts(2) + 1;
        case {'AU','au','Au','aU'}
            cWWCounts(3) = cWWCounts(3) + 1;
        case {'UA','ua','Ua','uA'}
            cWWCounts(4) = cWWCounts(4) + 1;
        case {'GU','gu','Gu','gU'}
            cWWCounts(5) = cWWCounts(5) + 1;
        case {'UG','ug','Ug','uG'}
            cWWCounts(6) = cWWCounts(6) + 1;
    end
    if ~isempty(breakPoint),
        cWW2 = strcat(sequence(breakPoint-1),sequence(breakPoint+1));
        switch cWW2,
            case {'CG','cg','Cg','cG'}
                cWWCounts(1) = cWWCounts(1) + 1;
            case {'GC','gc','Gc','gC'}
                cWWCounts(2) = cWWCounts(2) + 1;
            case {'AU','au','Au','aU'}
                cWWCounts(3) = cWWCounts(3) + 1;
            case {'UA','ua','Ua','uA'}
                cWWCounts(4) = cWWCounts(4) + 1;
            case {'GU','gu','Gu','gU'}
                cWWCounts(5) = cWWCounts(5) + 1;
            case {'UG','ug','Ug','uG'}
                cWWCounts(6) = cWWCounts(6) + 1;
        end
    else
        breakPoint = seqLength;
    end
    if breakPoint > 3,
        index0 = nuc2num(sequence(1));
        index = nuc2num(sequence(2));
        initialTransitionCounts(index0,index) = initialTransitionCounts(index0,index) +1;
        initialCounts(index) = initialCounts(index) + 1;
        coreCounts(index) = coreCounts(index) + 1;
    end
    if seqLength - breakPoint > 3,
        index0 = nuc2num(sequence(breakPoint+1));
        index = nuc2num(sequence(breakPoint+2));
        initialTransitionCounts(index0,index) = initialTransitionCounts(index0,index) +1;
        initialCounts(index) = initialCounts(index) + 1;
        coreCounts(index) = coreCounts(index) + 1;
    end
    if breakPoint > 4,
        index1 = nuc2num(sequence(2));
        for j = 3:breakPoint-2
            index2 = nuc2num(sequence(j));
            coreCounts(index2) = coreCounts(index2) + 1;
            transitionCounts(index1,index2) = transitionCounts(index1,index2) + 1;
            index1 = index2;
        end
    end
    if seqLength - breakPoint > 4,
        index1 = nuc2num(sequence(breakPoint+2));
        for j = breakPoint+3:seqLength-2
            index2 = nuc2num(sequence(j));
            coreCounts(index2) = coreCounts(index2) + 1;
            transitionCounts(index1,index2) = transitionCounts(index1,index2) + 1;
            index1 = index2;
        end
    end
end
end

function index = nuc2num(nuc)
index = -1;
switch nuc,
    case {'A','a'}
        index = 1;
    case {'C','c'}
        index = 2;
    case {'G','g'}
        index = 3;
    case {'U','u'}
        index = 4;
end
end