% pEditDistanceAll calculates the pairwise edit distances between the sequences in F and G,
% which are assumed to be of type loopType.
% It returns a length(F) x length(G) matrix of pairwise edit distances.
% location can be 'full' or 'core'
% It does not do any rotations.
% It returns one distance matrix.

function [D] = pEditDistanceAll(F,G,location)

    if isempty(G)
        G = F;
        self = 1;
    else
        self = 0;
    end

    if nargin < 3 || strcmp(location,'full')
        full = 1;
    else
        full = 0;
    end

    m = length(F);
    n = length(G);

    D = Inf * ones(m,n);

    if full == 1
        % loop over G first because that usually has more sequences
        for b = 1:n
            % split G sequence according to * character
            G_parts = strsplit(G(b).Sequence,'*');
            for a = 1:m
                % split F sequence according to * character
                F_parts = strsplit(F(a).Sequence,'*');
                d = 0;
                for k = 1:length(G_parts)
                    d = d + EditDist(F_parts{k},G_parts{k});
                end
                D(a,b) = d;

                % if d < 4
                %     fprintf('F sequence %20s G sequence %20s full distance %4d\n',F(a).Sequence,G(b).Sequence,d);
                % end
            end
        end
    else
        % loop over G first because that usually has more sequences
        for b = 1:n
            % split G sequence according to * character
            G_parts = strsplit(G(b).Sequence,'*');
            for k = 1:length(G_parts)
                gp = G_parts{k};
                % remove flanking bases
                G_parts{k} = gp(2:(end-1));
            end
            for a = 1:m
                % split F sequence according to * character
                F_parts = strsplit(F(a).Sequence,'*');
                d = 0;
                for k = 1:length(G_parts)
                    fp = F_parts{k};
                    % remove flanking bases
                    F_parts{k} = fp(2:(end-1));
                    d = d + EditDist(F_parts{k},G_parts{k});
                end
                D(a,b) = d;

                % if d < 4
                %     fs = F_parts{1};
                %     gs = G_parts{1};
                %     for k = 2:length(G_parts)
                %         fs = [fs '*' F_parts{k}];
                %         gs = [gs '*' G_parts{k}];
                %     end
                %     fprintf('F sequence %20s becomes %16s G sequence %20s becomes %16s core distance %4d\n',F(a).Sequence,fs,G(b).Sequence,gs,d);
                % end
            end
        end
    end
end