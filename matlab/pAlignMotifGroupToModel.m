% pAlignGroupToModel tells what each column in a motif group, and what each
% location for possible insertions between models, corresponds to in the model

% load 'MotifLibrary\IL_018.mat'
% JAR3D_path
% [filePaths,Node] = pWriteSingleJAR3DModel(Search,'IL_018',1)
% T = pAlignMotifGroupToModel(Search,Node,'IL_018','IL_018'); T'

function [T] = pAlignMotifGroupToModel(Search,Node,ModelName,MatFileName)

if nargin < 4,
  MatFileName = Search.Name;
end

MatFileName = strrep(MatFileName,'.mat','');

clear T
r = 1;

for n = 1:length(Node),
  if strcmp(Node(n).type,'Initial') && n > 1,
    if n < length(Node),
      if Node(n+1).LeftLetter == '*',
        i = min(Node(n).LeftIndex);
        T{r} = sprintf('%s_Column_%d-*_Insertion corresponds_to_JAR3D %s_Node_%d_Position_%d_Insertion', MatFileName, i-1, ModelName, n, 1);
        r = r + 1;
        i = Node(n).RightIndex;
        T{r} = sprintf('%s_Column_*-%d_Insertion corresponds_to_JAR3D %s_Node_%d_Position_%d_Insertion', MatFileName, i+1, ModelName, n, 2);
        r = r + 1;
      else
        if length(Node(n).leftLengthDist) > 1,  % has a non-zero insertion probability
          i = min(Node(n).LeftIndex);
          T{r} = sprintf('%s_Column_%d-%d_Insertion corresponds_to_JAR3D %s_Node_%d_Position_%d_Insertion', MatFileName, i-1, i, ModelName, n, 1);
          r = r + 1;
        end
        if length(Node(n).rightLengthDist) > 1, % has a non-zero insertion probability
          i = Node(n).RightIndex;
          T{r} = sprintf('%s_Column_%d-%d_Insertion corresponds_to_JAR3D %s_Node_%d_Position_%d_Insertion', MatFileName, i, i+1, ModelName, n, 2);
          r = r + 1;
        end
      end
    end
  end
  if strcmp(Node(n).type,'Fixed'),
    if length(Node(n).leftLengthDist) > 1,
      i = Node(n).LeftIndex;
      T{r} = sprintf('%s_Column_%d corresponds_to_JAR3D %s_Node_%d_Position_%d', MatFileName, i, ModelName, n, 1);
      r = r + 1;
    elseif length(Node(n).rightLengthDist) > 1,
      i = Node(n).RightIndex;
      T{r} = sprintf('%s_Column_%d corresponds_to_JAR3D %s_Node_%d_Position_%d', MatFileName, i, ModelName, n, 2);
      r = r + 1;
    end
  end
  if strcmp(Node(n).type,'Basepair'),
    i = Node(n).LeftIndex;
    T{r} = sprintf('%s_Column_%d corresponds_to_JAR3D %s_Node_%d_Position_%d', MatFileName, i, ModelName, n, 1);
    r = r + 1;
    if Node(n+1).LeftLetter == '*',
      T{r} = sprintf('%s_Column_%d-*_Insertion corresponds_to_JAR3D %s_Node_%d_Position_%d_Insertion', MatFileName, i, ModelName, n, 1);
      r = r + 1;
    else
      T{r} = sprintf('%s_Column_%d-%d_Insertion corresponds_to_JAR3D %s_Node_%d_Position_%d_Insertion', MatFileName, i, i+1, ModelName, n, 1);
      r = r + 1;
    end
    i = Node(n).RightIndex;
    T{r} = sprintf('%s_Column_%d corresponds_to_JAR3D %s_Node_%d_Position_%d', MatFileName, i, ModelName, n, 2);
    r = r + 1;
    if Node(n+1).LeftLetter == '*',
      T{r} = sprintf('%s_Column_*-%d_Insertion corresponds_to_JAR3D %s_Node_%d_Position_%d_Insertion', MatFileName, i, ModelName, n, 2);
      r = r + 1;
    else
      T{r} = sprintf('%s_Column_%d-%d_Insertion corresponds_to_JAR3D %s_Node_%d_Position_%d_Insertion', MatFileName, i-1, i, ModelName, n, 2);
      r = r + 1;
    end
  end
  if strcmp(Node(n).type,'Cluster')
    j = [Node(n).LeftIndex(Node(n).Left) Node(n).RightIndex(Node(n).Right)];
    for jj = 1:length(j),
      i = j(jj);
      T{r} = sprintf('%s_Column_%d corresponds_to_JAR3D %s_Node_%d_Position_%d', MatFileName, i, ModelName, n, jj);
      r = r + 1;
      if jj ~= length(Node(n).Left) && jj < length(j),
        T{r} = sprintf('%s_Column_%d-%d_Insertion corresponds_to_JAR3D %s_Node_%d_Position_%d-%d_Insertion', MatFileName, i, i+1, ModelName, n, jj, jj+1);
        r = r + 1;
      end
    end
  end
  if strcmp(Node(n).type,'Hairpin') && ~strcmp(Node(n).LeftLetter,'*'),
    j = Node(n).MiddleIndex;
    for jj = 1:length(j),
      i = j(jj);
      T{r} = sprintf('%s_Column_%d corresponds_to_JAR3D %s_Node_%d_Position_%d', MatFileName, i, ModelName, n, jj);
      r = r + 1;
      if jj < length(j),
        T{r} = sprintf('%s_Column_%d-%d_Insertion corresponds_to_JAR3D %s_Node_%d_Position_%d-%d_Insertion', MatFileName, i, i+1, ModelName, n, jj, jj+1);
        r = r + 1;
      end
    end
  end
end

