% pSmoothDistribution(values) returns a normalized distribution over the numbers in values, but that has more positive entries than the values would otherwise indicate

function [distn] = pSmoothDistribution(values)

    distn = zeros(1,max(values)+3);                             % make space to store probabilities
    L = length(values);

    for c = 1:L,
        if values(c) > 0,
            distn(values(c)+0) = distn(values(c)+0) + 1/(20*L); % one less than the value observed
        end
        if values(c) > -1,
          distn(values(c)+1) = distn(values(c)+1) + 1.00;       % the actual value observed
        end
        distn(values(c)+2) = distn(values(c)+2) + 1/(20*L);     % one more than the value observed
        distn(values(c)+3) = distn(values(c)+3) + 1/(400*L);    % two more than the value observed
    end

    w = find(distn);
    for i = min(w):max(w),
        if distn(i) == 0,
            distn(i) = 1/(8000*L);                              % fill in between observed values
        end
    end

    if all(values == 0),
        switch length(values),
        case 1
            distn = [0.99 0.01 0.0001];
        case 2
            distn = [0.992 0.008 0.00005];
        case 3
            distn = [0.994 0.006 0.00003];
        case 4
            distn = [0.996 0.004 0.00001];
        otherwise
            distn = [0.998 0.002 0.00001];
        end
    end

    distn = distn/sum(distn);

if 0 > 1,
    v = [0]
    pSmoothDistribution(v)
    v = [0 0]
    pSmoothDistribution(v)
    v = [0 0 0]
    pSmoothDistribution(v)
    v = [0 0 0 0]
    pSmoothDistribution(v)
    v = [0 0 0 0 0]
    pSmoothDistribution(v)
    v = [0 0 0 0 0 0]
    pSmoothDistribution(v)
    v = [1]
    pSmoothDistribution(v)
    v = [1 1]
    pSmoothDistribution(v)
    v = [1 1 1]
    pSmoothDistribution(v)
    v = [1 1 1 1]
    pSmoothDistribution(v)
    v = [1 1 1 1 1]
    pSmoothDistribution(v)
    v = [1 1 1 1 1 1]
    pSmoothDistribution(v)
    v = [0 1]
    pSmoothDistribution(v)
    v = [0 0 1 1]
    pSmoothDistribution(v)
    v = [0 0 0 1 1 1]
    pSmoothDistribution(v)
    v = [0 0 0 0 1 1 1 1]
    pSmoothDistribution(v)
    v = [0 1 1]
    pSmoothDistribution(v)
    v = [0 1 1 1]
    pSmoothDistribution(v)
end