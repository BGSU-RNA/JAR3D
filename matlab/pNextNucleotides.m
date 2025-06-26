% find next nucleotides on the strand ahead of a and behind B

function [aaa,BBB] = pNextNucleotides(a,B,Truncate,N,cdepth)

    if pStrandBreaksBetween(a,B,Truncate) == 0
        % hairpin
        aaa = min([N a+cdepth floor((a+B)/2)]);
        BBB = max([1 B-cdepth floor((a+B)/2)+1]);
    else
        % IL or junction
        aaa = min([N a+cdepth pNextTruncation(a,Truncate,N)]);
        BBB = max([1 B-cdepth pNextTruncation(B,Truncate,1)]);
    end
