% pColumnsForModel lists out what column each part of a JAR3D model should go into when displayed

% This program needs to call itself recursively to deal with junctions!!!

% load 'MotifLibrary\IL_018.mat'
% JAR3D_path
% [filePaths,Node] = pWriteSingleJAR3DModel(Search,'IL_018',1)
% T = pColumnsForModel(Search,Node,'TempModel','IL_018'); T'
% Sample output:
% IL_018_Node_3_Position_1 appears_in_column 10

function [T] = pColumnsForModelRecursive(Node,N)

clear T                                % build up lines from the left
r = 1;
clear U                                % build up lines from the right,
                                       % will be read backwards
s = 1;

n = N;
continu = 1;

while continu > 0,
  if strcmp(Node(n).type,'Initial')
    if length(Node(n).leftLengthDist) > 1,
      T{r} = sprintf('Node_%d_Position_1_Insertion', n);
      r = r + 1;
    end
    if length(Node(n).rightLengthDist) > 1,
      U{s} = sprintf('Node_%d_Position_2_Insertion', n);
      s = s + 1;
    end
  elseif strcmp(Node(n).type,'Fixed'),
    if length(Node(n).leftLengthDist) > 1,
      T{r} = sprintf('Node_%d_Position_1', n);
      r = r + 1;
    else
      U{s} = sprintf('Node_%d_Position_2', n);
      s = s + 1;
    end
  elseif strcmp(Node(n).type,'Basepair'),
    T{r} = sprintf('Node_%d_Position_1', n);
    r = r + 1;
    T{r} = sprintf('Node_%d_Position_1_Insertion', n);
    r = r + 1;
    U{s} = sprintf('Node_%d_Position_2', n);
    s = s + 1;
    U{s} = sprintf('Node_%d_Position_2_Insertion', n);
    s = s + 1;
  elseif strcmp(Node(n).type,'Cluster')
    k = 1;
    for j = 1:(length(Node(n).LeftIndex)-1),
      T{r} = sprintf('Node_%d_Position_%d', n, k);
      r = r + 1;
      T{r} = sprintf('Node_%d_Position_%d-%d_Insertion', n, k, k+1);
      r = r + 1;
      k = k + 1;
    end
    T{r} = sprintf('Node_%d_Position_%d', n, k);
    r = r + 1;
    k = k + 1;

    kk = length(Node(n).LeftIndex) + length(Node(n).RightIndex);

    for j = 1:(length(Node(n).RightIndex)-1),
      U{s} = sprintf('Node_%d_Position_%d', n, kk);
      s = s + 1;
      U{s} = sprintf('Node_%d_Position_%d-%d_Insertion', n, kk-1, kk);
      s = s + 1;
      kk = kk - 1;
    end
    U{s} = sprintf('Node_%d_Position_%d', n, k);
    s = s + 1;
  elseif strcmp(Node(n).type,'Hairpin')
    k = 1;
    for j = 1:(length(Node(n).MiddleIndex)-1),
      T{r} = sprintf('Node_%d_Position_%d', n, k);
      r = r + 1;
      T{r} = sprintf('Node_%d_Position_%d-%d_Insertion', n, k, k+1);
      r = r + 1;
      k = k + 1;
    end
    T{r} = sprintf('Node_%d_Position_%d', n, k);
    r = r + 1;

    T = [T fliplr(U)];                  % append, reverse entries of U
    continu = 0;                        % stop after a hairpin
  elseif strcmp(Node(n).type,'Junction')
    disp('Junction nodes have not been tested in pColumnsForModelRecursive');
    T = [T pColumnsForModelRecursive(Node,ModelName,Node(n).nextnode(1))];
    T = [T pColumnsForModelRecursive(Node,ModelName,Node(n).nextnode(2))];
    T = [T fliplr(U)];
    continu = 0;
  elseif strcmp(Node(n).type,'Alternative')
    disp('Alternative nodes are not set up in pColumnsForModelRecursive');
    continu = 0;
  end

  n = n + 1;
  if n > length(Node),
    continu = 0;
  end
end

