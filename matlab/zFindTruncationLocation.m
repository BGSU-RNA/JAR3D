% zFindTruncationLocation locates the location(s) of strand breaks

function [Truncate,index] = zFindTruncationLocation(loopType,Edge)

if strcmp(loopType,'IL'),
    index = find(diag(fix(abs(Edge)),1)==1);
    Truncate = index+1;
else
    index = [];
    Truncate = [];
end

