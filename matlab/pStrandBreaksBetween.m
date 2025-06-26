% pStrandBreaksBetween(a,B,Truncate) counts the number of strand breaks after a
% and before or at B

function [count] = pStrandBreaksBetween(a,B,Truncate)
    count = sum(B >= Truncate) - sum(a > Truncate);
end