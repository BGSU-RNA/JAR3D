% pInsertionDistribution finds the length and letter distribution following column a in Search.Candidates

function [LengthDist,LetterDist,insnumber,letter] = pInsertionDistribution(Search,a,UseIndex,Prior,Normalize);

	if nargin < 5,
		Normalize = 1;
	end

	if nargin < 4,
		Prior = [1 1 1 1];
	end

	if nargin < 3,
		UseIndex = length(Search.Candidates(:,1));
	end

	[L,N] = size(Search.Candidates);
	N = N - 1;
	L = length(UseIndex);

    insnumber = zeros(1,L);
    letter = zeros(size(Prior));

    for c = 1:L,
        insnumber(c) = double(abs(Search.Candidates(UseIndex(c),a+1) - Search.Candidates(UseIndex(c),a))) - 1;
        if insnumber(c) > 0,
            f = Search.Candidates(UseIndex(c),N+1);           % file number
            d = Search.Candidates(UseIndex(c),a);             % index within file
            for i = 1:insnumber(c),
                insb = Search.File(f).NT(d+i).Code;           % A=1, C=2, G=3, U=4
                if ~isempty(insb)
                    letter(insb) = letter(insb) + 1;
                end
            end
        end
    end

    LengthDist = pSmoothDistribution(insnumber);
    if Normalize == 1,
        LetterDist = (letter + Prior) / sum(letter + Prior);  % normalize
    else
    	LetterDist = [1 1 1 1];
    end
