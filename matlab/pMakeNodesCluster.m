% pMakeNodesCluster sets up cluster node of an appropriate size
% Later: identify clusters made up entirely of bases on the left
% or bases on the right.  Now, it combines them, even if they could
% be left distinct.  Build from the left and from the right, and
% when they overlap, coalesce them.

% Problem:  When the cluster is merely a set of interactions on the left
% strand or on the right strand, it still consumes one base on the other
% strand, even though that might not make sense.

% Note:  There are three numbering schemes:
% 1. Index from original file, so numbers could be very large.  a, b, aa, bb are on this scheme.
% 2. Numbering along each strand, from 1 to N on each strand.  ss, tt, zs, zt, xs, yt, RightNotInter are on this scheme.
% 3. Numbering the interacting bases in the cluster from 1 to M, regardless of strand.  rightnum, e, IBases are on this scheme.

if exist('zxsx')
  disp('Cluster *****************************************************');
end

%[a b B cdepth]
b = B;                                 % current base on right

amax = min(a+cdepth,floor((a+b)/2));   % how far to look on left
bmin = max(b-cdepth,floor((a+b)/2)+1); % how far to look on right

z = find((a < Truncate) .* (Truncate < amax));
if ~isempty(z),
  amax = Truncate(z)-1;                % don't look beyond truncation point
end

z = find((bmin < Truncate) .* (Truncate < b));
if ~isempty(z),
  bmin = Truncate(z)-1;                % don't look beyond truncation point
end

X = full(triu(G(a:amax,a:amax)));      % interactions on left strand
Y = full(triu(G(bmin:b,bmin:b)));      % interactions on right strand
Z = full(G(a:amax,bmin:b));            % interactions between left + right
[s,t] = size(Z);                       % X is s x s, Y is t x t

% ------------------------------------ determine extent of cluster interactions

ssa = max(find(X(1,:)));               % depth from a on left
ssb = max(find(Z(:,t)));               % nearest interaction to t 
if isempty(ssa), ssa = 1; end
if isempty(ssb), ssb = 1; end
ss = max(ssa, ssb);

tta = min(find(Z(1,:)));               % depth from b on left
ttb = min(find(Y(:,t)));               % depth from b on right
if isempty(tta), tta = t; end
if isempty(ttb), ttb = t; end
tt = min(tta, ttb);

while sum(sum(Z(1:ss,1:(tt-1)))) > 0 || ...
      sum(sum(Z((ss+1):s,tt:t))) > 0 || ...
      sum(sum(X(1:ss,(ss+1):s))) > 0 || ...
      sum(sum(Y(tt:t,1:(tt-1)))) > 0,
  ssa = max(find(sum(X(1:ss,:),1)));
  ssb = max(find(sum(Z(:,tt:t),2)));
  if isempty(ssa), ssa = 1; end
  if isempty(ssb), ssb = 1; end
  ss = max(ssa, ssb);

  tta = min(find(sum(Z(1:ss,:),1)));
  ttb = min(find(sum(Y(:,tt:t),2)));
  if isempty(tta), tta = t; end
  if isempty(ttb), ttb = t; end
  tt = min(tta, ttb);
end

aa = a - 1 + ss;                      % left extent of cluster
bb = b - (t - tt);                    % right extent of cluster

% ------------------------------------- Set up basics of cluster node

n=n+1;
Node(n).Delete       = 0.001;            % deletion probability
Node(n).type         = 'Cluster';        % node type
Node(n).nextnode     = n+1;              % index of next node in tree
Node(n).LeftIndex    = [a:aa];           % full range of indices
Node(n).LeftLetter   = cat(2,File.NT(a:aa).Base);
Node(n).RightLetter  = cat(2,File.NT(bb:b).Base);
Node(n).RightIndex   = [bb:b];

AllIndices = [a:aa bb:b];

Node(n).Comment = [' // Cluster node ' File.NT(Node(n).LeftIndex(1)).Base File.NT(Node(n).LeftIndex(1)).Number ':' File.NT(Node(n).LeftIndex(end)).Base, File.NT(Node(n).LeftIndex(end)).Number ' and ' File.NT(Node(n).RightIndex(1)).Base File.NT(Node(n).RightIndex(1)).Number ':' File.NT(Node(n).RightIndex(end)).Base File.NT(Node(n).RightIndex(end)).Number];

if Verbose > 0,
  fprintf('%3d Cluster   %s:%s %s:%s\n', n, File.NT(a).Number, File.NT(aa).Number, File.NT(bb).Number, File.NT(b).Number);
end

% -------------------- identify which bases are interacting on left and right

X = X(1:ss,1:ss);          % reduce interactions to this cluster only
Y = Y(tt:t,tt:t);
Z = Z(1:ss,tt:t);

zs = find(sum(Z,2)');      % bases on left interacting with right
zt = find(sum(Z,1));       % bases on right interacting with left
xs = find(sum(X+X',1));    % left with left
yt = find(sum(Y+Y',1));    % right with right

Node(n).LeftNotInter = setdiff(1:length(Node(n).LeftIndex),union(zs,xs));
Node(n).RightNotInter = setdiff(1:length(Node(n).RightIndex),union(zt,yt));

if insertionconserved == 1,  % treat non-interacting bases as part of the cluster
  zs = 1:ss;                 % all are thought of as interacting
  zt = 1:length(Y(1,:));
  xs = 1:ss;
  yt = 1:length(Y(1,:));
end

leftinter = union(zs,xs);            % list of interacting on left
leftnum   = [];
leftnum(leftinter) = 1:length(leftinter); 
                          % sequential numbering of interacting bases

rightinter = union(zt,yt);
rightnum   = [];
rightnum(rightinter) = 1:length(rightinter);
                          % sequential numbering of interacting bases
rightnum(find(rightnum)) = rightnum(find(rightnum))+length(leftinter);  % shift to numbering for all bases in cluster

% --------------- list insertion possibilities and corresponding probabilities

zsxs = union(zs,xs);              % all interacting bases on the left
if isempty(zsxs),                 % happens when no inter across
  zsxs = 1;
end

ztyt = union(zt,yt);              % all interacting bases on the right
if isempty(ztyt),
  ztyt = 1;
end

Node(n).Left(1,:)  = zsxs;                  % which bases interact
Node(n).Right(1,:) = ztyt;                  % which bases interact

Node(n).PIns = [0.00001 0.99999];           % when no previous state
Node(n).P    = ones(17,1) * Node(n).PIns;   % transition probabilities - old

% ------------- create a list of insertions according to what is observed here

e = [Node(n).Left(1,:) Node(n).Left(1,end)+Node(n).Right(1,:)];
                                            % indices of bases used, left to right

d = diff(e);                                % diffs in positions used
h = find(d>1);                              % where insertions occur
cc = 1;                                     % counter for insertions
for aaa = 1:length(h),                      % loop through insertions
  Node(n).Insertion(cc).Position = h(aaa);

  ID = pMakeNodesInsertionDist(cat(2,File.NT(AllIndices(h(aaa)+1):AllIndices(h(aaa)+d(h(aaa))-1)).Code),'Cluster',Normalize);

  if d(aaa) == 2,                           % exactly one insertion
disp('Checking for interactions made by single inserted base');
    R = pAdjustSubsProb(File,AllIndices(h(aaa)+1),[],ID.LetterDist,method);
  end

  Node(n).Insertion(cc).LengthDist = ID.LengthDist;
  Node(n).Insertion(cc).LetterDist = ID.LetterDist;

  Node(n).InsertionComment{aaa} = [' // Insertion between ' File.NT(AllIndices(h(aaa))).Base File.NT(AllIndices(h(aaa))).Number ' and ' File.NT(AllIndices(h(aaa)+d(h(aaa)))).Base File.NT(AllIndices(h(aaa)+d(h(aaa)))).Number];

  if Verbose > 0,
    fprintf('    %d insertions (%s) between %s%s and %s%s\n', ...
    d(h(aaa))-1, ...
    cat(2,File.NT((AllIndices(h(aaa))+1):(AllIndices(h(aaa)+d(h(aaa)))-1)).Base), ...
                 File.NT(AllIndices(h(aaa))).Base, ...
                 File.NT(AllIndices(h(aaa))).Number, ...
                 File.NT(AllIndices(h(aaa)+d(h(aaa)))).Base, ...
                 File.NT(AllIndices(h(aaa)+d(h(aaa)))).Number);
  end
  cc = cc + 1;
end

% -------------------------------------- create a list of interacting bases

% ---- interactions between left and right
K = 1;                                % counter
[i,j] = find(Z);                      % interacting pairs
for k = 1:length(i),                  % loop through them
    Node(n).IBases(K,:) = [leftnum(i(k)) rightnum(j(k))];
    i1 = i(k) + a - 1;                % index of first base
    i2 = j(k) + bb - 1;               % index of second
    Node(n).InterIndices(K,:) = [i1 i2];

    Node(n).InteractionComment{K} = [ ' // Cluster Interaction ' File.NT(i1).Base File.NT(i1).Number ' - ' File.NT(i2).Base File.NT(i2).Number ' ' zEdgeText(File.Edge(i1,i2))];

    if Verbose > 0,
      fprintf('    LR Inter  %4s %4s %c%c %s\n', File.NT(i1).Number, File.NT(i2).Number, File.NT(i1).Base, File.NT(i2).Base, zEdgeText(File.Edge(i1,i2)));
%              fprintf('%d %d\n', Node(n).IBases(K,1), Node(n).IBases(K,2));
    end

    Node(n).SubsProb(:,:,K) = pIsoScore(File.Edge(i1,i2), ...
       File.NT(i1).Code, File.NT(i2).Code,method,ExemplarIDI,ExemplarFreq,Normalize);

    Node(1).Edge(i1,i2) = File.Edge(i1,i2);
    K  = K + 1;
end

% ---- interactions between left and left
[i,j] = find(X);                     % interacting pairs
for k = 1:length(i),                 % loop through them
    Node(n).IBases(K,:) = [leftnum(i(k)) leftnum(j(k))];
    i1 = i(k) + a - 1;             % index of first base
    i2 = j(k) + a - 1;             % index of second
    Node(n).InterIndices(K,:) = [i1 i2];

    Node(n).InteractionComment{K} = [ ' // Cluster Interaction ' File.NT(i1).Base File.NT(i1).Number ' - ' File.NT(i2).Base File.NT(i2).Number ' ' zEdgeText(File.Edge(i1,i2))];

    if Verbose > 0,
      fprintf('    LL Inter  %4s %4s %c%c %s\n', File.NT(i1).Number, File.NT(i2).Number, File.NT(i1).Base, File.NT(i2).Base, zEdgeText(File.Edge(i1,i2)));
%fprintf('%d %d\n', Node(n).IBases(K,1), Node(n).IBases(K,2));
    end

    Node(n).SubsProb(:,:,K) = pIsoScore(File.Edge(i1,i2), ...
 File.NT(i1).Code, File.NT(i2).Code,method,ExemplarIDI,ExemplarFreq,Normalize);

    Node(1).Edge(i1,i2) = File.Edge(i1,i2);
    K  = K + 1;
end

% ---- interactions between right and right
[i,j] = find(Y);                      % interacting pairs
for k = 1:length(i),                  % loop through them
    Node(n).IBases(K,:) = [rightnum(i(k)) rightnum(j(k))];
    i1 = i(k) + bb - 1;               % index of first base
    i2 = j(k) + bb - 1;               % index of second
    Node(n).InterIndices(K,:) = [i1 i2];

    Node(n).InteractionComment{K} = [ ' // Cluster Interaction ' File.NT(i1).Base File.NT(i1).Number ' - ' File.NT(i2).Base File.NT(i2).Number ' ' zEdgeText(File.Edge(i1,i2))];

    if Verbose > 0,
      fprintf('    RR Inter  %4s %4s %c%c %s\n', File.NT(i1).Number, File.NT(i2).Number, File.NT(i1).Base, File.NT(i2).Base, zEdgeText(File.Edge(i1,i2)));
%fprintf('%d %d\n', Node(n).IBases(K,1), Node(n).IBases(K,2));
    end

    Node(n).SubsProb(:,:,K) = pIsoScore(File.Edge(i1,i2), ...
     File.NT(i1).Code, File.NT(i2).Code,method,ExemplarIDI,ExemplarFreq,Normalize);
    Node(1).Edge(i1,i2) = File.Edge(i1,i2);    % record interaction
    K  = K + 1;
end

% ---- inserted bases which are treated as conserved insertions and modeled as interacting with themselves

for i = Node(n).LeftNotInter,
  K = length(Node(n).IBases(:,1))+1;
  Node(n).IBases(K,:) = [i i];
  i1 = i + a - 1;
  Node(n).InterIndices(K,:) = [i1 i1];
  Node(n).InteractionComment{K} = [' // Left strand conserved insertion ' File.NT(i1).Base File.NT(i1).Number];
  M = eye(4);
  M(File.NT(i1).Code,File.NT(i1).Code) = 4;
  M = M / sum(sum(M));
  Node(n).SubsProb(:,:,K) = M;            % favor the observed base
  
  if Verbose > 0,
    fprintf('    Left strand conserved insertion %c%4s at position %d\n', File.NT(i1).Base, File.NT(i1).Number, i);
  end
end

for i = Node(n).RightNotInter,
  K = length(Node(n).IBases(:,1))+1;

  Node(n).IBases(K,:) = [i i] + length(leftinter); % shift from right strand numbering to whole motif numbering
  i1 = i + bb - 1;                        % index in original file
  Node(n).InterIndices(K,:) = [i1 i1];
  Node(n).InteractionComment{K} = [' // Right strand conserved insertion ' File.NT(i1).Base File.NT(i1).Number];
  M = eye(4);
  M(File.NT(i1).Code,File.NT(i1).Code) = 4;
  M = M / sum(sum(M));
  Node(n).SubsProb(:,:,K) = M;            % favor the observed base

  if Verbose > 0,
    fprintf('    Right strand conserved insertion %c%4s at right strand position %d\n', File.NT(i1).Base, File.NT(i1).Number, i);
  end

end

% ----------------------------------------- Adjust subs probs for LR, BPh inter

if AdjustSubsForLR == 1 && (TertiaryFreeNode == 0 || Extension < 2),
  for K = 1:length(Node(n).IBases(:,1)),    % run through basepairs
    i1 = Node(n).InterIndices(K,1);         % left base
    i2 = Node(n).InterIndices(K,2);         % right base
    P  = Node(n).SubsProb(:,:,K);           % current substitution probs
    if Verbose > 1,
      NT1 = File.NT(i1);
      NT2 = File.NT(i2);
      fprintf('    Possibly adjusting substitution probabilities for %4s %4s %s%s %s\n', NT1.Number, NT2.Number, NT1.Base, NT2.Base, zEdgeText(File.Edge(i1,i2)));
    end
    Node(n).SubsProb(:,:,K) = pAdjustSubsProb(File,i1,i2,P,method,AllIndices);
  end
end

a = aa + 1;                           % current base on left
                                      % skip over rest of cluster
B = bb - 1;                           % current base on right

a = min(a,B);                         % just in case!

% calculate normalization constant
%if Normalize == 1,
    Node(n).NormCons = pClusterNorm(Node(n).InterIndices,Node(n).SubsProb,Node(n).LeftIndex,Node(n).RightIndex);
%else
%    Node(n).NormCons = pClusterNorm(Node(n).InterIndices,Node(n).SubsProb,Node(n).LeftIndex,Node(n).RightIndex);
%    Node(n).NormCons = .66;             % pClusterNorm only works for normalized matrices, temp fix!
%end
    
%         psum = 0;
%         numBases = length(Node(n).Left) + length(Node(n).Right);
% 		if numBases < 11,
%             % System.out.println("ClusterNode.Normalize "+numBases);
%             for i = 1:4^numBases,
%                 for j = 1:numBases
%                     li = floor((i-1)/(4^(numBases-j)));
%                     li = mod(li,4)+1;
%                     code(j) = li;
%                 end
%                 [numInter,dum] = size(Node(n).InterIndices);
%                 Cols(1:length(Node(n).Left)) = Node(n).LeftIndex(Node(n).Left);
%                 Cols(length(Node(n).Left)+1:length(Node(n).Left)+length(Node(n).Right)) = Node(n).RightIndex(Node(n).Right);
%                 % score codes[] according to the various interactions
%                 prob = 1;
%                 for j = 1:numInter,
%                     i1 = find(Cols == Node(n).InterIndices(j,1));
%                     i2 = find(Cols == Node(n).InterIndices(j,2));
%                     prob = prob*Node(n).SubsProb(code(i1),code(i2),j);
%                 end
%                 psum = psum+prob;     % running total of probabilities
%                 %			System.out.println(z);
%                 % push a 1 into codes and carry
%             end
%             Node(n).NormCons = psum;
%         else
%             Node(n).NormCons = 1;
%         end
%		System.out.println("ClusterNode.Normalize normalized a node: " + z);
