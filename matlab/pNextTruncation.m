% pNextTruncation(a,Truncate,N) finds the last value in the same strand as a is in,
% on the way to N

function [c] = pNextTruncation(a,Truncate,N)

    % make sure that Truncate is a row vector so the function returns just one value
    if size(Truncate,1) > 1
        Truncate = Truncate';
    end

    c = N;
    if N == 1
        % look earlier in the sequence, toward 1
        for t = Truncate
            if t <= a
                c = t;
            end
        end
    else
        % look later in the sequence, toward N
        for t = Truncate
            if t > a
                c = t-1;  % last index before the next strand starts
                break;
            end
        end
    end
end
