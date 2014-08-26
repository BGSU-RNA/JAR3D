% pConservedInteractionList(Search,SubsProb) returns a cell array of
% conserved pairwise interactions.  When SubsProb = 1, it appends a
% list of substitution probabilities.

function [Text] = pConservedInteractionList(Search,SubsProb,Verbose)

if nargin < 3,
    Verbose = 0;
end

[i,j,k] = find(triu(Search.Edge));        % use consensus interactions

[a,b] = sortrows([i j],[1 2]);

i = i(b);
j = j(b);
k = k(b);

Lett = 'ACGU';

r = 1;
Text{1} = '';                              % initialize, so something gets passed back

for s = 1:length(b),                       % loop through conserved interactions
    aa = i(s);
    bb = j(s);
    if isfield(Search,'SubsProb')
        if isempty(Search.SubsProb{aa,bb}),
            M = ones(4,4)/16;
            if abs(k(s)) < 20,
                if Verbose >= 0,
                    fprintf('pConservedInteractionList did not find a substitution probability matrix for %s between %d and %d\n', zEdgeText(k(s)), aa, bb);
                end
            end
        else
            M = Search.SubsProb{aa,bb};
        end
    else
        M = ones(4,4)/16;
        if abs(k(s)) < 20,
            if Verbose >= 0,
                fprintf('pConservedInteractionList did not find a substitution probability matrix for %s between %d and %d\n', zEdgeText(k(s)), aa, bb);
            end
        end
    end

    sp = '';
    if SubsProb > 0,
        for x = 1:4,
            for y = 1:4,
                sp = [sp sprintf('%c%c %0.4f ', Lett(x), Lett(y), M(x,y))];
            end
        end
    end

    if SubsProb == 0 || abs(k(s)) < 20,
        Text{r} = sprintf('%d %d %s', i(s), j(s), zEdgeText(k(s)));
        Text{r} = [Text{r} sp];
        r = r + 1;
        if Verbose > 0,
            fprintf('pConservedInteractionList: %s\n', Text{r});
        end
    end
end

if SubsProb == 0,

    [i,j,k] = find(Search.BPh);        % use consensus interactions

    for s = 1:length(i),                       % loop through conserved interactions
        aa = i(s);
        bb = j(s);
        Text{r} = sprintf('%d %d %s', i(s), j(s), strrep(zBasePhosphateText(k(s)),' ',''));
        Text{r} = [Text{r} sp];
        r = r + 1;
    end

    [i,j,k] = find(Search.BR);        % use consensus interactions

    for s = 1:length(i),                       % loop through conserved interactions
        aa = i(s);
        bb = j(s);
        Text{r} = sprintf('%d %d %s', i(s), j(s), strrep(zBaseRiboseText(k(s)),' ',''));
        Text{r} = [Text{r} sp];
        r = r + 1;
    end

end