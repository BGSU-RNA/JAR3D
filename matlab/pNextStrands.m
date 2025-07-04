% find next nucleotides on the strand after a and before B
% for junction situations only

function [NSL,NSR] = pNextStrands(a,B,Truncate,N)

    r = pNextTruncation(a,Truncate,N);  % end of a strand
    s = r+1;  % start of next strand after a
    t = pNextTruncation(s,Truncate,N);

    NSL = s:t;

    z = pNextTruncation(B,Truncate,1);  % beginning of B strand
    y = z-1;                            % end of strand before B
    x = pNextTruncation(y,Truncate,1);  % start of strand before B

    NSR = x:y;
end