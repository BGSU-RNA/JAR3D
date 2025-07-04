function [cWWCounts, initialCounts, transitionCounts, coreCounts, initialTransitionCounts] = pExtractSequenceCounts(Text)

    % Text is a cell array alternating between header and sequence

    cWWCounts = zeros(6,1);           %CG,GC,AU,UA,GU,UG
    initialCounts = zeros(4,1);       %A,C,G,U
    coreCounts = zeros(4,1);          %A,C,G,U
    transitionCounts = zeros(4,4);    %(A,C,G,U)^2
    initialTransitionCounts = zeros(4,4);    %(A,C,G,U)^2
    numSeq = floor(length(Text)/2);
    for i = 1:numSeq
        sequence = Text{2*i};
        seqLength = length(sequence);
        breakPoint = strfind(sequence, '*');  % locations of * characters

        % outer cWW pair, which is present in every loop
        cWW1 = strcat(sequence(1),sequence(seqLength));
        switch upper(cWW1)
            case 'CG'
                cWWCounts(1) = cWWCounts(1) + 1;
            case 'GC'
                cWWCounts(2) = cWWCounts(2) + 1;
            case 'AU'
                cWWCounts(3) = cWWCounts(3) + 1;
            case 'UA'
                cWWCounts(4) = cWWCounts(4) + 1;
            case 'GU'
                cWWCounts(5) = cWWCounts(5) + 1;
            case 'UG'
                cWWCounts(6) = cWWCounts(6) + 1;
        end

        for bp = breakPoint
            % count cWW pairs across the * characters for IL, J3, etc.
            cWW2 = strcat(sequence(bp-1),sequence(bp+1));
            switch upper(cWW2)
                case 'CG'
                    cWWCounts(1) = cWWCounts(1) + 1;
                case 'GC'
                    cWWCounts(2) = cWWCounts(2) + 1;
                case 'AU'
                    cWWCounts(3) = cWWCounts(3) + 1;
                case 'UA'
                    cWWCounts(4) = cWWCounts(4) + 1;
                case 'GU'
                    cWWCounts(5) = cWWCounts(5) + 1;
                case 'UG'
                    cWWCounts(6) = cWWCounts(6) + 1;
            end
        end

        starts = [1 breakPoint+1];   % starting point of every strand
        ends = [breakPoint-1 seqLength];   % ending points of every strand

        % loop over strands
        for j = 1:length(starts)
            seq = sequence(starts(j):ends(j));

            if length(seq) > 3
                % process start of seq
                index0 = nuc2num(seq(1));
                index = nuc2num(seq(2));
                initialTransitionCounts(index0,index) = initialTransitionCounts(index0,index) + 1;
                initialCounts(index) = initialCounts(index) + 1;
                coreCounts(index) = coreCounts(index) + 1;
            end

            if length(seq) > 4
                % process transitions within interior of seq
                index1 = nuc2num(seq(2));
                for k = 3:length(seq)-1
                    index2 = nuc2num(seq(k));
                    coreCounts(index2) = coreCounts(index2) + 1;
                    transitionCounts(index1,index2) = transitionCounts(index1,index2) + 1;
                    index1 = index2;
                end
            end
        end

        if 0 > 1
            % earlier code that handled HL and IL
            % special cases depending on how long the strands are ...
            % previously breakpoint was just a single number, * location for IL, end of strand for HL
            % looks like we need to just process each strand separately using this code
            % did the second strand of IL ever get processed?
            if breakPoint > 3
                % start of first strand if it is long enough
                index0 = nuc2num(sequence(1));
                index = nuc2num(sequence(2));
                initialTransitionCounts(index0,index) = initialTransitionCounts(index0,index) +1;
                initialCounts(index) = initialCounts(index) + 1;
                coreCounts(index) = coreCounts(index) + 1;
            end
            if seqLength - breakPoint > 3
                % start of second strand if it is long enough
                index0 = nuc2num(sequence(breakPoint+1));
                index = nuc2num(sequence(breakPoint+2));
                initialTransitionCounts(index0,index) = initialTransitionCounts(index0,index) +1;
                initialCounts(index) = initialCounts(index) + 1;
                coreCounts(index) = coreCounts(index) + 1;
            end
            if breakPoint > 4
                % transitions in core of first strand if it is long enough
                index1 = nuc2num(sequence(2));
                for j = 3:breakPoint-2
                    index2 = nuc2num(sequence(j));
                    coreCounts(index2) = coreCounts(index2) + 1;
                    transitionCounts(index1,index2) = transitionCounts(index1,index2) + 1;
                    index1 = index2;
                end
            end
            if seqLength - breakPoint > 4
                % transitions in core of second strand if it is long enough
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
end

function index = nuc2num(nuc)
    index = -1;
    switch upper(nuc)
        case 'A'
            index = 1;
        case 'C'
            index = 2;
        case 'G'
            index = 3;
        case 'U'
            index = 4;
        otherwise
            fprintf('pExtractSequenceCounts: Unknown nucleotide %s\n',nuc);
    end
end