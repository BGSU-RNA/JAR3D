% pConsensusPairSubstitution(a,b,f,File,F,L,Search) looks at the letter pairs corresponding to nucleotides a and b of the query motif in Search, uses interactions in F as the consensus, 

% The argument Noncanonical being set to 1 means that this pair should set low probabilities for CG, GC, AU, and UA base combinations, because this pair is the first non-canonical pair in a motif.

function [Score] = pConsensusPairSubstitution(a,b,f,File,F,Search,Param,Normalize,Noncanonical,Prior)

if nargin < 9,
  Noncanonical = [];
end

Verbose = Param(1);

Score = zeros(4,4);                       % ready to sum IsoScores for this pair 
Count = zeros(4,4);                       % tally observed basepairs here

[L,N] = size(Search.Candidates);          % L = num instances; N = num NT
N = N - 1;                                % number of nucleotides

if a ~= b,                                % two different bases

  method = 2;                             % default method for assigning pair subst probs

  if length(Param) > 1,
    method  = Param(2);
  end

  load PairExemplars

  cl = fix(F.Edge(a,b));                  % consensus pair or stack between a and b

  if any(cl == [-1 -2 -7 -8 -12 -22 -23 -101 -102 -107 -108 -112 -122 -123]),
    cl = -cl;                             % symmetric pair, use positive class
  end

  if Verbose > 0,
    fprintf('pConsensusPairSubstitution: Consensus is %4s\n', zEdgeText(cl));
  end

  % --------------------------------- Accumulate 4x4 matrices from instances

  for c = 1:L,                            % loop through candidates
    i   = Search.Candidates(c,a);         % index 1 of pair in candidate
    j   = Search.Candidates(c,b);         % index 2 of pair in candidate
    NT1 = File(f(c)).NT(i);               % retrieve the first nucleotide
    NT2 = File(f(c)).NT(j);               % retrieve the second nucleotide

    g = File(f(c)).Edge(i,j);             % actual interaction for this pair

    if Verbose > 0,
      fprintf('pConsensusPairSubstitution: File %4s has %s%4s and %s%4s making %4s\n', File(f(c)).Filename, NT1.Base, NT1.Number, NT2.Base, NT2.Number, zEdgeText(g));
    end

    if ismember([a b],Noncanonical,'rows'),
%      fprintf('pConsensusPairSubstitution: File %4s has %s%4s and %s%4s making %4s\n', File(f(c)).Filename, NT1.Base, NT1.Number, NT2.Base, NT2.Number, zEdgeText(g));
%      pause
    end

    g = fix(g);                           % actual interaction

    if any(g == [-1 -2 -7 -8 -12 -22 -23 -101 -102 -107 -108 -112 -122 -123]),
      g = -g;                             % symmetric pair, use positive
    end

    if g == cl || g == sign(cl)*(100+abs(cl)), % near or true match to consensus
      newScore = pIsoScore(cl,NT1.Base,NT2.Base,method,ExemplarIDI,ExemplarFreq,Normalize);
                                               % use consensus edge


      v = Search.Candidates(c,1:N);            % indices of nucleotides in current instance

      for e = 1:3,                             % loop through BPh and BR edges
        w = (Search.BPh(a,:) == e) .* File(f(c)).BasePhosphate(i,v);
%        w(a) = 0;                              % exclude self interaction
        if Param(10) > 0 && sum(w) > 0,        % scoring with BPh interactions
 if w(a) > 0,
   fprintf('pConsensusPairSubstitution:  New 4x4 substitution matrix for instance %d of %d, pair %d,%d for nucleotide %d (%s) making self BPh with edge %d\n', c, L, a,b,a, NT1.Base,e);
 end
          newScore(NT1.Code,:) = Param(10)*newScore(NT1.Code,:);  % conservation of this base
          newScore = newScore/sum(sum(newScore)); % normalize
        elseif Param(11) > 0,                  % scoring with BR interactions
          w = (Search.BR(a,:) == e) .* File(f(c)).BaseRibose(i,v);
          w(a) = 0;                            % exclude self interaction
          w(b) = 0;                            % exclude BR when a and b are making a basepair
          if sum(w) > 0,                       % a makes conserved BPh and one in this instance
            k = find(w);
            k = k(1);
%  fprintf('pConsensusPairSubstitution:  New 4x4 substitution matrix for instance %d of %d, pair %d,%d for nucleotide %d (%s) making BR with nucleotide %d using edge %d\n', c, L, a,b,a, NT1.Base,k,e);
            newScore(NT1.Code,:) = Param(11)*newScore(NT1.Code,:);  % conservation of this base
            newScore = newScore/sum(sum(newScore)); % normalize
          end
        end

        w = (Search.BPh(b,:) == e) .* File(f(c)).BasePhosphate(j,v);
%        w(b) = 0;                              % exclude self interaction
        if Param(10) > 0 && sum(w) > 0,        % b makes conserved BPh and one in this instance
 if w(b) > 0,
%  fprintf('pConsensusPairSubstitution:  New 4x4 substitution matrix for instance %d of %d, pair %d,%d for nucleotide %d (%s) making BPh with edge %d\n', c, L, a,b,b, NT2.Base,e);
 end
          newScore(:,NT2.Code) = Param(10)*newScore(:,NT2.Code);  % conservation of this base
          newScore = newScore/sum(sum(newScore)); % normalize
        elseif Param(11) > 0,                  % scoring with BR interactions
          w = (Search.BR(b,:) == e) .* File(f(c)).BaseRibose(j,v);
          w(b) = 0;                            % exclude self interaction
          w(a) = 0;                            % exclude BR when a and b are making a basepair
          if sum(w) > 0,                       % b makes conserved BPh and one in this instance
            k = find(w);
            k = k(1);
%  fprintf('pConsensusPairSubstitution:  New 4x4 substitution matrix for instance %d of %d, pair %d,%d for nucleotide %d (%s) making BR with nucleotide %d using edge %d\n', c, L, a,b,b, NT2.Base,k,e);
            newScore(:,NT2.Code) = Param(11)*newScore(:,NT2.Code);  % 50% conservation of this base
            newScore = newScore/sum(sum(newScore)); % normalize
          end
        end
      end

      Score = Score + newScore;
      if Verbose > 0,
        fprintf(', using it with isostericity\n');
      end
    else                                     % does not match consensus
      Score(NT1.Code,NT2.Code) = Score(NT1.Code,NT2.Code) + 0.1;  % count this, but less
      if Verbose > 0,
        fprintf(', using it without isostericity\n');
      end
    end

    Count(NT1.Code,NT2.Code) = Count(NT1.Code,NT2.Code) + 1;
  end

  Score = Score / sum(sum(Score));           % normalize
  Count = Count / sum(sum(Count));           % normalize

  if Verbose > 0,
    Score
    Count
  end

  Score = (1-(L/(L+Param(9)))) * Score + (L/(L+Param(9))) * Count;   % weighted average of scoring methods

else                                      % conserved but non-basepairing position

  Count = zeros(size(Prior));                     % place to store actual counts

  for c = 1:L,                            % loop through candidates
    i   = Search.Candidates(c,a);         % index 1 of pair in candidate
    NT1 = File(f(c)).NT(i);               % retrieve the first nucleotide

    if Verbose > 0,
      fprintf('pConsensusPairSubstitution: File %4s has non-basepairing core base %s%4s\n', File(f(c)).Filename, NT1.Base, NT1.Number);
    end

    Count(1,NT1.Code) = Count(1,NT1.Code) + 1;    
  end

  Pr = Prior;

  for e = 1:3,                                        % loop through edges
      if Param(10) > 0 && sum(Search.BPh(a,:)==e) > 0,
          Pr = Pr / Param(10);                  % weaken the prior distribution
      elseif Param(11) > 0 && sum(Search.BR(a,:)==e) > 0,
          Pr = Pr / Param(11);                  % weaken the prior distribution
      end
  end

  Count = Count + Pr;                       % add the prior
  Score = diag(Count / sum(Count));         % normalize

end
