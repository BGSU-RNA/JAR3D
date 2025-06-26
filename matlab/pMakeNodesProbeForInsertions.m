% probe for insertions on left and right strands, set insertion probabilities
% if current node is a cluster, always insert an initial node after that
% there should be no interaction between a and B

aa = a;                                  % initial location on left strand
BB = B;                                  % initial location on right strand

movea = 1;                               % plan to advance a until we see a reason not to
moveB = 1;

if pStrandBreaksBetween(a,B,Truncate) > 1
  % there is a junction ahead, don't skip ahead super greedily
  % this section is for junction motifs, not for entire molecules

  % indices of last nucleotides on current left and right strands
  aaa = pNextTruncation(a,Truncate,N);
  BBB = pNextTruncation(B,Truncate,1);

  % nucleotides on next strand(s) (NSL and NSR are the same for J3)
  [NSL,NSR] = pNextStrands(a,B,Truncate,N);

  while movea == 1
    if a == aaa
      % reached the end of this strand
      movea = 0;
    elseif sum(sum(H(a,NSL))) > sum(sum(H(a,BBB:B)))
      % more interactions with NSL than right strand, stop
      movea = 0;
    elseif sum(sum(H(a,BBB:B))) == 0 && sum(sum(H(a,(a+1):aaa))) == 0
      % a makes no interactions
      fprintf('    Insertion on left  %s\n', File.NT(a).ID);
      a = a + 1;
    else
      movea = 0;
    end
  end

  while moveB == 1
    if B == BBB
      % reached the end of this strand
      moveB = 0;
    elseif sum(sum(H(B,NSR))) > sum(sum(H(B,a:aaa)))
      % more interactions with NSL than right strand, stop
      moveB = 0;
    elseif sum(sum(H(B,a:aaa))) == 0 && sum(sum(H(B,BBB:(B-1)))) == 0
      % B makes no interactions
      fprintf('    Insertion on right %s\n', File.NT(B).ID);
      B = B - 1;
    else
      moveB = 0;
    end
  end

else
  % HL and IL probing
  % ---------------------------------------- Carefully written probing algorithm!
  %                                          Appropriate for HL and IL
  [aaa,BBB] = pNextNucleotides(a,B,Truncate,N,cdepth);

  % ---------------------------------------- Easy insertions on the left
  while sum(sum(H(a,a:B))) == 0 && a < aaa   % a makes no basepairs w/in loop
    if Verbose > 0
      fprintf('    Insertion on left  %s\n', File.NT(a).ID);
    end
    a = a + 1;
  end

  % ---------------------------------------- Easy insertions on the right
  while sum(sum(H(a:B,B))) == 0 && B > BBB   % B makes no basepairs w/in loop
    if Verbose > 0
      fprintf('    Insertion on right %s\n', File.NT(B).ID);
    end
    B = B - 1;                               % next base on right
  end

  % ---------------------------------------- a interacts, but skip anyway?
  LS  = (a+1):aaa;                           % left strand not including a
  RS  = BBB:B;                               % right strand including B
  aB  = a:B;                                 % indices from a to B
  na  = find(H(a,aB));                       % nested interaction(s) a makes
  nB  = find(H(B,aB));                       % nested interaction(s) B makes

  if ~isempty(na) && Verbose > 0
    fprintf('pMakeNodesProbeForInsertions:  a %s interacts with %s\n', File.NT(a).ID, File.NT(aB(max(na))).ID);
  end

  if ~isempty(nB) && Verbose > 0
    fprintf('pMakeNodesProbeForInsertions:  B %s interacts with %s\n', File.NT(B).ID, File.NT(aB(min(nB))).ID);
  end

  if a == aaa
    % reached the end of the strand
    movea = 0;
  elseif sum(sum(H((a+1):aB(max(na)),(a+1):(aB(max(na)))))) > 0
    % after a and up to the last interaction partner of a, there are basepairs
    % that will always happen in an IL if a makes a basepair with the other strand
    movea = 0;
    if Verbose > 0
      disp('pMakeNodesProbeForInsertions:  Found that the interaction a makes contains a stem');
    end
  elseif ~isempty(nB)                         % B makes a nested interaction
    if min(abs(nonzeros(G(a,aB(na))))) < min(abs(nonzeros(G(B,aB(nB)))))
      movea = 0;                              % a makes a more important one
      if Verbose > 0
        disp('pMakeNodesProbeForInsertions:  Both a and B make nested interactions, but a makes a more important one');
      end
    else
      movea = 1;                              % B makes a more important one
      if Verbose > 0
        disp('pMakeNodesProbeForInsertions:  Both a and B make nested interactions, but B makes a more important one');
      end
    end
  elseif any(abs(G(a,aB(na)))==1)
    % a makes a cWW with something up to B
    movea = 0;
  end

  % accumulate interactions on the left
  % could just set movea = 1 and check the conditions at the start of the while loop
  while sum(abs(G(a,[LS RS]))) == 0 && movea == 1
                                          % no interaction w/in this loop
    if Verbose > 0
      fprintf('    Insertion on left$ %s\n', File.NT(a).ID);
    end
    a = a + 1;                              % next base on left
    [aaa,BBB] = pNextNucleotides(a,B,Truncate,N,cdepth);
    LS  = (a+1):aaa;                        % left strand
    RS  = BBB:B;                            % right strand
    aB  = a:B;                              % indices from a to B
    na  = find(H(a,aB));                       % nested interaction(s) a makes
    nB  = find(H(B,aB));                       % nested interaction(s) B makes

    if a == aaa
      % reached the end of the strand
      movea = 0;
    elseif sum(sum(H((a+1):aB(max(na)),(a+1):(aB(max(na)))))) > 0
      movea = 0;
    elseif ~isempty(nB)                        % B makes a nested interaction
      if min(abs(nonzeros(G(a,aB(na))))) < min(abs(nonzeros(G(B,aB(nB)))))
        movea = 0;
      end
    elseif any(abs(G(a,aB(na)))==1)
      movea = 0;
    end
  end

  % ---------------------------------------- More easy insertions on the right

  while sum(sum(G(a:B,B))) == 0 && B > BBB   % B makes no interactions w/in loop
    if Verbose > 0
      fprintf('    Insertion on right %s\n', File.NT(B).ID);
    end
    B = B - 1;                               % next base on right
  end

  % --------------------------------------- B interacts, but skip anyway?

  [aaa,BBB] = pNextNucleotides(a,B,Truncate,N,cdepth);
  LS  = a:aaa;                               % left strand including a
  RS  = BBB:(B-1);                           % right strand not including B
  aB  = a:B;                                 % indices from a to B
  na  = find(H(a,aB));                       % nested interaction(s) a makes
  nB  = find(H(B,aB));                       % nested interaction(s) B makes

  if ~isempty(nB) && Verbose > 0
    fprintf('pMakeNodesProbeForInsertions:  B %s interacts with %s\n', File.NT(B).ID, File.NT(aB(max(nB))).ID);
  end

  if B == BBB
    % reached the end of the strand
    moveB = 0;
  elseif sum(sum(H(aB(min(nB)):(B-1),aB(min(nB)):(B-1)))) > 0
    moveB = 0;
    if Verbose > 0
      disp('pMakeNodesProbeForInsertions:  Found that the interaction B makes contains a stem');
    end
  elseif ~isempty(nB)                        % B makes a nested interaction
    if min(abs(nonzeros(G(B,aB(nB))))) < min(abs(nonzeros(G(a,aB(na))))),
      moveB = 0;
      if Verbose > 0
        disp('pMakeNodesProbeForInsertions:  Both a and B make nested interactions, but B makes a more important one');
      end
    end
  elseif any(abs(G(B,aB(nB)))==1)
    moveB = 0;
  end

  while sum((G([LS RS],B))) == 0 && moveB == 1
                                          % no interaction w/in this loop
    if Verbose > 0
      fprintf('    Insertion      %4s %s\n', File.NT(B).Number, File.NT(B).Base);
    end
    B = B - 1;                               % next base on right
    [aaa,BBB] = pNextNucleotides(a,B,Truncate,N,cdepth);
    LS  = a:aaa;                             % left strand
    RS  = BBB:(B-1);                         % right strand
    aB  = a:B;                            % indices from a to B
    na  = find(H(a,aB));                       % nested interaction(s) a makes
    nB  = find(H(B,aB));                       % nested interaction(s) B makes

    if B == BBB
      % reached the end of the strand
      moveB = 0;
    elseif sum(sum(H(aB(min(nB)):(B-1),aB(min(nB)):(B-1)))) > 0
      moveB = 0;
    elseif ~isempty(nB)                        % B makes a nested interaction
      if min(abs(nonzeros(G(B,aB(nB))))) < min(abs(nonzeros(G(a,aB(na)))))
        moveB = 0;
      end
    elseif any(abs(G(B,aB(nB)))==1)
      moveB = 0;
    end
  end
end

LeftIns  = a-aa;                            % number of insertions on the left
RightIns = BB-B;                            % number of insertions on the right

% ---------------------------------------------- Set insertion parameters
% a and B have been moved as much as possible, now set parameters in the right place

uniform = [1 1 1 1]'/4;

if strcmp(Node(n).type,'Basepair')    % add insertions to basepair
  i1 = Node(n).LeftIndex+1;           % index of first inserted base on left
  LeftCodes = cat(2,File.NT(i1:(a-1)).Code);
  ID = pMakeNodesInsertionDist(LeftCodes,'Basepair',Normalize);
  if length(LeftCodes) == 1 && AdjustSubsForLR == 1  % exactly one insertion
    R = pAdjustSubsProb(File,i1,[],uniform,method); % adjust for LR basepairs
    ID.LetterDist = R;
  end

  Node(n).leftLengthDist  = ID.LengthDist;
  Node(n).leftLetterDist  = ID.LetterDist;

  if length(LeftCodes) > 0 && Verbose > 1
    fprintf('pMakeNodesProbeForInsertions:  LeftIns %d  length(LeftCodes) %d\n', LeftIns, length(LeftCodes));
  end

  i1 = Node(n).RightIndex-1;           % index of first inserted base on right
  RightCodes = cat(2,File.NT((B+1):(Node(n).RightIndex-1)).Code);
  ID = pMakeNodesInsertionDist(RightCodes,'Basepair',Normalize);
  if length(RightCodes) == 1 && AdjustSubsForLR == 1 % exactly one insertion
    R = pAdjustSubsProb(File,i1,[],uniform,method); % adjust for LR basepairs
    ID.LetterDist = R;
    disp('Adjusted for LR interactions');
  end
  Node(n).rightLengthDist = ID.LengthDist;
  Node(n).rightLetterDist = ID.LetterDist;

  if length(RightCodes) > 0 && Verbose > 1
    fprintf('pMakeNodesProbeForInsertions:  RightIns %d  length(RightCodes) %d\n', RightIns, length(RightCodes));
  end

  Node(n).LeftLetter  = [Node(n).LeftLetter cat(2,File.NT((Node(n).LeftIndex+1):(a-1)).Base)];
  Node(n).RightLetter = [cat(2,File.NT((B+1):(Node(n).RightIndex-1)).Base) Node(n).RightLetter];

% ----------------------------------------------- Insertions for Initial node

elseif strcmp(Node(n).type,'Initial')
  i1 = Node(n).LeftIndex;               % index of first inserted base on left
  LeftCodes = cat(2,File.NT(i1:(a-1)).Code);
  ID = pMakeNodesInsertionDist(LeftCodes,'Initial',Normalize);
  if length(LeftCodes) == 1 && AdjustSubsForLR == 1,  % exactly one insertion
    R = pAdjustSubsProb(File,i1,[],uniform,method);% adjust for LR basepairs
    ID.LetterDist = R;
  end
  Node(n).leftLengthDist  = ID.LengthDist;
  Node(n).leftLetterDist  = ID.LetterDist;

  i1 = Node(n).RightIndex;           % index of first inserted base on right
  RightCodes = cat(2,File.NT((B+1):i1).Code);
  ID = pMakeNodesInsertionDist(RightCodes,'Initial',Normalize);
  if length(RightCodes) == 1 && AdjustSubsForLR == 1 % exactly one insertion
    R = pAdjustSubsProb(File,i1,[],uniform,method); % adjust for LR basepairs
    ID.LetterDist = R;
  end
  Node(n).rightLengthDist = ID.LengthDist;
  Node(n).rightLetterDist = ID.LetterDist;

  Node(n).LeftLetter  = cat(2,File.NT(Node(n).LeftIndex:(a-1)).Base);
  Node(n).RightLetter = cat(2,File.NT((B+1):Node(n).RightIndex).Base);

% ----------------------------------------------- Insertion node after cluster

% elseif strcmp(Node(n).type,'Cluster') && (LeftIns > 0 || RightIns > 0)
elseif strcmp(Node(n).type,'Cluster') % always put an initial node after a cluster node
  n = n + 1;                          % initial node after cluster
  Node(n).type = 'Initial';
  Node(n).nextnode  = n+1;            % index of next node in tree
  Node(n).LeftIndex  = aa;
  Node(n).RightIndex = BB;
  Node(n).LeftLetter  = cat(2,File.NT(aa:(aa+LeftIns-1)).Base);
  Node(n).RightLetter = cat(2,File.NT((BB-RightIns+1):BB).Base);

  i1 = Node(n).LeftIndex;               % index of first inserted base on left
  LeftCodes = cat(2,File.NT(i1:(i1+LeftIns-1)).Code);
  ID = pMakeNodesInsertionDist(LeftCodes,'Initial',Normalize);
  if length(LeftCodes) == 1 && AdjustSubsForLR == 1  % exactly one insertion
    R = pAdjustSubsProb(File,i1,[],uniform,method);% adjust for LR basepairs
    ID.LetterDist = R;
  end
  Node(n).leftLengthDist  = ID.LengthDist;
  Node(n).leftLetterDist  = ID.LetterDist;

  i1 = Node(n).RightIndex;           % index of first inserted base on right
  RightCodes = cat(2,File.NT((BB-RightIns+1):i1).Code);
  ID = pMakeNodesInsertionDist(RightCodes,'Initial',Normalize);
  if length(RightCodes) == 1 && AdjustSubsForLR == 1 % exactly one insertion
    R = pAdjustSubsProb(File,i1,[],uniform,method); % adjust for LR basepairs
    ID.LetterDist = R;
  end
  Node(n).rightLengthDist = ID.LengthDist;
  Node(n).rightLetterDist = ID.LetterDist;

end

% ---------------------------------------------- Add comments for initial node

if strcmp(Node(n).type,'Initial')
  Node(n).Comment = [' // Initial node ' File.NT(aa).ID ' - ' File.NT(BB).ID ' from Node ' Node(n).id];

  if Verbose > 0
    fprintf('%3d Initial   %4s (%d insertion) and %4s (%d insertion)\n', n, File.NT(Node(n).LeftIndex).Number, LeftIns, File.NT(Node(n).RightIndex).Number, RightIns);
  end
end


