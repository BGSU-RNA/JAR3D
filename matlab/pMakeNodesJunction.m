% pMakeNodesJunction decides how to split up strands around a junction
% and then makes a Node that splits off the first new branch

% current position in the molecule is a on the left and B on the right
r = a;                               % last index in the current strand before the junction
rr= a;                               % start of current loop

% this approach works for junction motifs, but not for junctions in molecules
aaa = pNextTruncation(a,Truncate,N);
BBB = pNextTruncation(B,Truncate,1);

% next strand left, next strand right
[NSL,NSR] = pNextStrands(a,B,Truncate,N);

% find last nt on next strand left that interacts with current left strand
s = NSL(end);
while sum(sum(G(s,a:aaa))) == 0 && s > NSL(1)
  s = s - 1;
end

t = min([s+1 NSL(end) B]);           % next after that, keep it in range
u = B;                               % end of current known loop

if Verbose > 0
    fprintf('pMakeNodesJunction: a=%s r=%s s=%s t=%s u=B=%s\n', File.NT(a).ID, File.NT(r).ID, File.NT(s).ID, File.NT(t).ID, File.NT(B).ID);
    fprintf('pMakeNodesJunction: a:s count %d\n', full(sum(sum(CanonicalcWW(r:s,r:s)==1))));
    fprintf('pMakeNodesJunction: t:B count %d\n', full(sum(sum(CanonicalcWW(t:u,t:u)==1))));
end

% if length(Truncate) > 1 && (sum(sum(CanonicalcWW(r:s,r:s)==1)) > 0) && (sum(sum(CanonicalcWW(t:u,t:u)==1)) > 0)
% not modeling an IL or HL and ...
% there are nested canonical cWW pairs between r and s and between t and u
% use cWW pairs to avoid pairs between i and i+1 making junctions

if Verbose > 0
    fprintf('pMakeNodesJunction: Found nested interactions between %s and %s and between %s and %s\n', File.NT(r).ID, File.NT(s).ID, File.NT(t).ID, File.NT(u).ID);
    fprintf('pMakeNodesJunction: about to make a junction: a=%s s=%s t=%s B=%s\n', File.NT(a).ID, File.NT(s).ID, File.NT(t).ID, File.NT(B).ID);
end

jcdepth = 5;   % how far to look for junction clusters ... which don't exist

fprintf('pMakeNodesJunction: first branch starts at   %s\n', File.NT(r).ID);
fprintf('pMakeNodesJunction: first branch ends around %s\n', File.NT(s).ID);
fprintf('pMakeNodesJunction: second branch starts at  %s\n', File.NT(t).ID);
fprintf('pMakeNodesJunction: junction ends at         %s\n', File.NT(u).ID);

rplus  = min([r+jcdepth s-1 pNextTruncation(r,Truncate,N)]);  % how far to look on the strand past r
sminus = max([s-jcdepth r+1 pNextTruncation(s,Truncate,1)]);  % how far to look on the strand before s
tplus  = min([t+jcdepth u-1 pNextTruncation(t,Truncate,N)]);  % how far to look on the strand past t
uminus = max([u-jcdepth t+1 pNextTruncation(u,Truncate,1)]);  % how far to look on the strand before u

fprintf('pMakeNodesJunction: first branch strand start goes up to  %s\n', File.NT(rplus).ID);
fprintf('pMakeNodesJunction: first branch strand end goes down to  %s\n', File.NT(sminus).ID);
fprintf('pMakeNodesJunction: second branch strand start goes up to %s\n', File.NT(tplus).ID);
fprintf('pMakeNodesJunction: last strand of junction goes down to  %s\n', File.NT(uminus).ID);

% count interactions beyond r and s that would cross the junction
C1 = full(sum(sum(H(r:rplus,t:tplus))));   % junction strands  1-3
C2 = full(sum(sum(H(r:rplus,uminus:u))));   % junction strands  1-4
C3 = full(sum(sum(H(sminus:s,t:tplus))));   % junction strands  2-3
C4 = full(sum(sum(H(sminus:s,uminus:u))));   % junction strands  2-4

if C1+C2+C3+C4 > 0 && Verbose > 0
  disp('pMakeNodesJunction: Junction includes some crossing pairs that will be lost. C1, C2, C3, C4:');
  disp([C1 C2 C3 C4])
end

% [i,j,k] = find(G(r:rplus,t:tplus));   % junction strands  1-3
% kk = k;
% [i,j,k] = find(G(r:rplus,uminus:u));   % junction strands  1-4
% kk = [kk; k];
% [i,j,k] = find(G(sminus:s,t:tplus));   % junction strands  2-3
% kk = [kk; k];
% [i,j,k] = find(G(sminus:s,uminus:u));   % junction strands  2-4
% kk = [kk; k];

% Node(1).JunctionDeletion = [Node(1).JunctionDeletion; kk];

% if Verbose > 1
%   disp('pMakeNodesJunction: This is how many basepairs are being lost to simplify the junction:');
%   disp(full([C1 C2 C3 C4]))
% end

% ----- Remove interactions that require junction cluster node -----
% Because those have not been developed yet
% Note:  this is only happening between the two loops identified so
% far, but there may be more loops, and between them these interactions
% are not being removed!  Maybe this is not actually a problem.
% A little study tells me that in E. coli 16S, no nested pairs
% were removed by this.

if 0 > 1
  if C1 > 0
    G(r:rplus,t:tplus) = 0*G(r:rplus,t:tplus);
  end
  if C2 > 0
    G(r:rplus,uminus:u) = 0*G(r:rplus,uminus:u);
  end
  if C3 > 0
    G(sminus:s,t:tplus) = 0*G(sminus:s,t:tplus);
  end
  if C4 > 0
    G(sminus:s,uminus:u) = 0*G(sminus:s,uminus:u);
  end

  C1 = full(sum(sum(G(r:rplus,t:tplus))));   % junction strands  1-3
  C2 = full(sum(sum(G(r:rplus,uminus:u))));   % junction strands  1-4
  C3 = full(sum(sum(G(sminus:s,t:tplus))));   % junction strands  2-3
  C4 = full(sum(sum(G(sminus:s,uminus:u))));   % junction strands  2-4
end

% probably we can just build the junction without zeroing out anything
C1 = 0;
C2 = 0;
C3 = 0;
C4 = 0;

if C1+C2+C3+C4 == 0        % no interaction across junction, so it's a plain junction
  % r to s will form one branch, t to u will form the others and they can be split further later

  junc = [];              % each row tells indices where branches start and end

  junc = [[r s]; [t u]];

  % note below, r+1 might equal s, then this would not be what you want
  % identify all of the branches coming off this junction,
  % each starts at r and ends at s, which keep changing as you go to the next branch

  % previously this was used to split up the junction
  % but now we may already have a good method
  branchesRemain = 0;
  while branchesRemain
  % while (r < pNextTruncation(r,Truncate,N)) && ...
  %       (s > pNextTruncation(s,Truncate,1)) && ...
  %       (sum(sum(H((r+1):(s-1),(r+1):(s-1)) > 0)) > 0) && ...
  %       (sum(sum(H((s+1):(u),(s+1):(u)) > 0)) > 0)

    % probe for start of next block of nested pairs,
    % s+1 to t-1 is a strand of the junction
    tStrandEnd = pNextTruncation(t,Truncate,N);
    while sum(H(t,(t+1):B) > 0) == 0 && t < u && t < tStrandEnd
      % if t does not make a nested pair later in this junction
      t = t + 1;
    end

    fprintf('pMakeNodesJunction: first branch could end as late as %s\n', File.NT(t-1).ID);

    tPartner = Interact{t}.Index(1);  % this had better be after t!
    if tPartner <= t
      fprintf('!!!!!!!!!!! tPartner is not after t, but %s\n', File.NT(tPartner).ID);
    end

    % decide how to split the strand between s and t between the branches

    % Note: rr is defined to equal a before pMakeNodesJunction is called!
    currentBranch = unique([rr:min(r+cdepth,s) max(rr,s-cdepth):s]); % indices of the current branch
    nextBranch = unique([t:min(t+cdepth,u) max(tPartner-cdepth,t):u]); % indices of the next branch

    % choose whichever endpoint has the larger number of interactions.
    % But it would be better yet to split the stand optimally, as with 556:567 in 2avy

    if sum(sum(G(currentBranch,(s+1):(t-1))>0)) >= sum(sum(G((s+1):(t-1),nextBranch))),
      % append start and end of the current branch to junc
      junc = [junc; [rr t-1]];      % current branch includes entire next strand
      r = t;
    else
      % append start and end of the current branch to junc
      junc = [junc; [rr s]];        % current branch ends at s
      r = s+1;
    end

    if Verbose > 0
      st = junc(end,1);         % start of this branch
      en = junc(end,2);         % end of this branch
      fprintf('Junction branch starts at %s ends at %s\n', File.NT(st).ID, File.NT(en).ID);
    end

    rr = r;                         % copy of r, for probing forward
    s  = Interact{t}.Index(1);      % the far side of the next branch
    t  = s + 1;                     % one position beyond s

    % I did not think about the next line, it is a gift from co-pilot
    if rr >= s || t > u || t > B || s > B || rr > B
      branchesRemain = 0;          % no more branches to find
    end
    junc = [junc; [r u]];             % append limits of this last branch

  end

  NL = length(junc(:,1));           % number of branches = number of rows appended

  id = fix(10000*rand);             % random id for large structures with many junctions

  if Verbose > 0
    fprintf('\nJunction with %d loops, call it J%d\n', NL,id);
    for ln = 1:NL
      fprintf('Actual Branch %d of junction J%d - Nucleotides %s to %s, length %3d\n',ln,id,File.NT(junc(ln,1)).ID,File.NT(junc(ln,2)).ID,junc(ln,2)+1-junc(ln,1));
    end
  end

  n = n+1;                              % move to next node
  Node(n).type       = 'Junction';      % junction node
  Node(n).LeftIndex  = a;
  Node(n).RightIndex = B;
  Node(n).NumLoops   = NL;              % not sure about this
  Node(n).id         = ['J' num2str(id)];
  Node(n).Comment = ['// Junction node ' File.NT(Node(n).LeftIndex).Base File.NT(Node(n).LeftIndex).Number ' - ' File.NT(Node(n).RightIndex).Base File.NT(Node(n).RightIndex).Number ' ID ' Node(n).id];

  if Verbose > 0
    fprintf('%3d Junction\n', n);
  end

  jn = n;                               % index of this node

  if NL == 2                            % exactly two branches,

    Node(jn).numbranches = 2;           % more reliable than tracking length(nextnode)
    for ln = 1:NL
      if Verbose > 0
        fprintf('\n');
        fprintf('2B Branch %d of %d of junction J%d - Nucleotides %5s to %5s, length %3d\n',ln,NL,id,File.NT(junc(ln,1)).Number,File.NT(junc(ln,2)).Number,junc(ln,2)+1-junc(ln,1));
      end
      nn = length(Node) + 1;
      Node(jn).nextnode(ln) =  length(Node)+1;
      Node = pMakeNodes(File,Param,junc(ln,1),junc(ln,2),Truncate,Data,Node,n);
      Node(nn).id         = ['J' num2str(id)];
      n = length(Node);
    end

  elseif NL > 2                           % more than two branches

    NN = ceil(NL/2);                      % # branches for 1st child

    if Verbose > 0
      fprintf('\n');
      fprintf('Branch %d of %d of junction J%d - Nucleotides %5s to %5s, length %3d\n',1,2,id,File.NT(junc(1,1)).Number,File.NT(junc(NN,2)).Number,junc(NN,2)+1-junc(1,1));
    end

    nn = length(Node) + 1;
    Node(jn).nextnode(1) = length(Node) + 1;
    Node = pMakeNodes(File,Param,junc(1,1),junc(NN,2),Truncate,Data,Node,n);
    Node(nn).id         = ['J' num2str(id)];
    n = length(Node);

    if Verbose > 0
      fprintf('\n');
      fprintf('Branch %d of %d of junction J%d - Nucleotides %5s to %5s, length %3d\n',2,2,id,File.NT(junc(NN+1,1)).Number,File.NT(junc(NL,2)).Number,junc(NL,2)+1-junc(NN+1,1));
    end

    % collect all remaining branches and process, maybe with a new junction
    % that way, each junction node only has two branches and life is simpler
    nn = length(Node) + 1;
    Node(jn).nextnode(2) = length(Node) + 1;
    Node = pMakeNodes(File,Param,junc(NN+1,1),junc(NL,2),Truncate,Data,Node,n);
    Node(nn).id         = ['J' num2str(id)];
    Node(jn).numbranches = length(Node(jn.nextnode));

  end

  % Node(n).Comment = ['// Junction node ' File.NT(Node(n).LeftIndex).Base File.NT(Node(n).LeftIndex).Number ' - ' File.NT(Node(n).RightIndex).Base File.NT(Node(n).RightIndex).Number ' ID ' Node(n).id];

else

  % find extent of interactions between loops

  t = b;

  rr = r + jcdepth;
  while (sum(G(rr,[t:tplus uminus:u])) == 0) && (rr > r),
    rr = rr - 1;
  end

  ss = s - jcdepth;
  while (sum(G(ss,[t:tplus uminus:u])) == 0) && (ss < s),
    ss = ss + 1;
  end

  tt = t + jcdepth;
  while (sum(G([r:rplus sminus:s],tt)) == 0) && (tt > t),
    tt = tt - 1;
  end

  uu = u - jcdepth;
  while (sum(G([r:rplus sminus:s],uu)) == 0) && (uu < u),
    uu = uu + 1;
  end

  % second, extent of additional interactions within loops

  currentBranch = r + jcdepth;
  while (sum(G(currentBranch,[r:rr ss:s])) == 0) && (currentBranch > rr),
    currentBranch = currentBranch - 1;
  end

  sss = s - jcdepth;
  while (sum(G(sss,[r:rr ss:s])) == 0) && (sss < ss),
    sss = sss + 1;
  end

  nextBranch = t + jcdepth;
  while (sum(G([t:tt uu:u],nextBranch)) == 0) && (nextBranch > tt),
    nextBranch = nextBranch - 1;
  end

  uuu = u - jcdepth;
  while (sum(G([t:tt uu:u],uuu)) == 0) && (uuu < uu),
    uuu = uuu + 1;
  end

  n = n + 1;
  Node(n).type        = 'JunctionCluster';  %
  Node(n).LeftIndex   = [r:currentBranch];
  Node(n).MiddleIndex = [sss:nextBranch];
  Node(n).RightIndex  = [uuu:u];

%        Node(n).Left(1,:)   = union(zs,xs);
%        Node(n).Middle(1,:) = fliplr(union(zt,yt)); % correct?
%        Node(n).Right(1,:)  = fliplr(union(zt,yt)); % correct?

  Node(n).LIP = [1];
  Node(n).MIP = [1];
  Node(n).RIP = [1];

  % add additional insertion combinations and probabilities here!
  % add scores for the various basepairs here!

  if Verbose > 0,
    fprintf('%3d Junction Cluster %4s %4s %4s %4s %4s %4s\n', n, File.NT(r).Number, File.NT(currentBranch).Number, File.NT(sss).Number, File.NT(nextBranch).Number, File.NT(uuu).Number, File.NT(u).Number);
    fprintf('================================================================================================================\n');
  end

  r = currentBranch;
  s = sss;
  t = nextBranch;
  u = uuu;

  Node(n).nextnode(1) =  n+1;          % index of next node in tree
  Node = pMakeNodes(File,Param,r,s,Truncate,Data,Node,n);

  Node(n).nextnode(2)  = length(Node)+1;
  Node = pMakeNodes(File,Param,t,u,Truncate,Data,Node,length(Node));
end                                  % junction cluster

Node(n).P    = [0.05*ones(17,1) 0.95*ones(17,1)]; % state to state transitions

Node(n).PIns = [0.05 0.95];   % when no previous state

EndLoop = 1;

