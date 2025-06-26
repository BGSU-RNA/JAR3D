% Add a hairpin node

n = n + 1;
Node(n).type = 'Hairpin';
MI                  = a:B;                % empty if a > B
Node(n).MiddleIndex = MI;                 % fixed bases in the hairpin
Node(n).LeftIndex   = a;           % left index of hairpin bases
Node(n).RightIndex  = B;           % right index of hairpin bases
Node(n).LeftLetter  = cat(2,File.NT(MI).Base);
Node(n).RightLetter = '';
Node(n).P           = ones(17,1);
Node(n).PIns        = 1;

if TertiaryFreeNode == 0 || Extension < 2   % this stem has long-range inter
                                            % or we are using regular hairpins
  Node(n).subtype     = 'XXXX';             % revise this later!

  % ----------------------------------------- model interactions
  [i,j] = find(triu(G(MI,MI)));
  for k = 1:length(i)
    Node(n).IBases(k,:) = [i(k) j(k)];      % these bases interact
    i1 = i(k) + a - 1;                      % index of first base
    i2 = j(k) + a - 1;                      % index of second base
    Node(n).SubsProb(:,:,k) = pIsoScore(File.Edge(i1,i2), ...
       File.NT(i1).Code, File.NT(i2).Code,method,ExemplarIDI,ExemplarFreq,Normalize);
    Node(n).InteractionComment{k} = [ ' // Hairpin Interaction ' File.NT(i1).Base File.NT(i1).Number ' - ' File.NT(i2).Base File.NT(i2).Number ' ' zEdgeText(File.Edge(i1,i2))];
  end

  NotBasepaired = setdiff(1:length(MI),union(i,j));

  for i = NotBasepaired,
    [K,KK] = size(Node(n).IBases);
    K = K + 1;

    Node(n).IBases(K,:) = [i i];
    i1 = i + a - 1;                        % index in original file
    Node(n).InterIndices(K,:) = [i1 i1];
    Node(n).InteractionComment{K} = [' // Hairpin conserved non-basepairing position ' File.NT(i1).Base File.NT(i1).Number];
    M = eye(4);
    M(File.NT(i1).Code,File.NT(i1).Code) = 4;
    Node(n).SubsProb(:,:,K) = M / sum(sum(M));            % favor the observed base

    if Verbose > 0
      fprintf('    Right strand conserved insertion %c%4s at right strand position %d\n', File.NT(i1).Base, File.NT(i1).Number, i);
    end

  end

  Node(n).Comment = [ ' // Hairpin node ' File.NT(a).Base File.NT(a).Number ':' File.NT(B).Base File.NT(B).Number];

  if Verbose > 0
    fprintf('%3d Hairpin   %s:%s %s\n', n, File.NT(a).Number, File.NT(B).Number, cat(2,File.NT(MI).Base));
  end

else                                        % no long-range interactions
  Node(n).subtype     = 'Vague';                  % revise this later!

  % The code below *assumes* that the hairpin has at least 3 nucleotides,
  % which it may not.  It is attempting to model a vague closing pair,
  % and presumably the previous basepair is forced to be cWW by some
  % other program

  % --------------------------------------- % interaction between 1 and 3

  Node(n).IBases      = [1 3];              % fixed bases 1 and 3 interact

  Node(n).InteractionComment{1} = [' // Hairpin on extensible stem ' File.NT(a).Number ' - ' File.NT(B).Number ' ' cat(2,File.NT(MI).Base)];

  P = [1 1 1 10; 1 1 10 1; 8 10 1 1; 10 4 8 1];   % typical closing pairs

  Node(n).SubsProb(:,:,1) = P / sum(sum(P));      % normalize subs probs

  if any(nonzeros(fix(File.Edge(MI,MI))) == -6) ...  % tSW, as in UNCG
     || any(nonzeros(fix(File.Edge(MI,MI))) == 10), ... %tSH as in GNRA
     Node(1).Edge(MI,MI) = File.Edge(MI,MI);       % as if we got this one
  end

  Node(n).Insertion(1).Position = 2;
  Node(n).Insertion(1).LengthDist = [0.4 0.4 0.15 0.05];
  if Normalize == 1,
    Node(n).Insertion(1).LetterDist = [1 1 1 1]/4;
  else
    Node(n).Insertion(1).LetterDist = [1 1 1 1];
  end
  Node(n).InsertionComment{1} = ' // Insertion in vague hairpin';

  Node(n).Comment = [ ' // Hairpin node on extensible stem ' File.NT(a).Base File.NT(a).Number ':' File.NT(B).Base File.NT(B).Number];

  if Verbose > 0,
    fprintf('%3d Hairpin   %s:%s %s  Extensible stem\n', n, File.NT(a).Number, File.NT(B).Number, cat(2,File.NT(MI).Base));
  end

end

EndLoop = 1;

% ----------------------------------------- Look for long-range interactions

if Verbose > 1,
  LR = File.Edge .* (File.Crossing > 2);
  [i,j,e] = find(LR(MI,:));
  i = i + a - 1;
  [yy,ii] = sort(abs(e));
  i = i(ii);
  j = j(ii);
  e = e(ii);

  if length(i) > 0
    fprintf('Hairpin has these long-range interactions: ===================\n');
    for z = 1:length(i)
      NT1 = File.NT(i(z));
      NT2 = File.NT(j(z));
      fprintf('Pair %s %s%5s_%s - %s%5s_%s %s\n', File.Filename, NT1.Base,NT1.Number,NT1.Chain,NT2.Base,NT2.Number,NT2.Chain, zEdgeText(File.Edge(i(z),j(z))));
    end
  elseif Verbose > 1
    fprintf('Hairpin has no long-range interactions ***********************\n');
  end

  c = cat(1,File.NT(1:length(File.NT)).Center); % nucleotide centers
  cent = mean(c);

  if Verbose > 1
    fprintf('Distance from center of molecule: %8.4f\n', norm(cent-File.NT(a).Center));
  end

  if Verbose > 2
    fprintf('Press any key to continue\n');
    pause
  end
end
%Calculate normalization constant
Node(n).InterIndices = Node(n).MiddleIndex(Node(n).IBases);
if ~isempty(Node(n).InterIndices)
    Left = Node(n).MiddleIndex(1:length(Node(n).MiddleIndex)-1);
    Node(n).NormCons = pClusterNorm(Node(n).InterIndices,Node(n).SubsProb,Left,Node(n).RightIndex);
else
    Node(n).NormCons = 1;
end