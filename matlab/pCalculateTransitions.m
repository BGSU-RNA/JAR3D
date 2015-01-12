ModelTestName = '';
DiagnosticTestName = '';

Prior = [.5 .5 .5 .5 0];              % Prior distribution for insertion bases

% ----------------------------------- initialize transistion matrix counts
cWWCounts = zeros(6,1);           %CG,GC,AU,UA,GU,UG
initialCounts = zeros(4,1);       %A,C,G,U
coreCounts = zeros(4,1);          %A,C,G,U
transitionCounts = zeros(4,4);    %(A,C,G,U)^2
initialTransitionCounts = zeros(4,4);    %(A,C,G,U)^2
TransitionFile = [OutputPath filesep 'transitions.mat'];

for m = 1:length(Filenames),
  MotifName = Filenames(m).name;
  load([MotifLibraryPath filesep MotifName '.mat']);
  [Search,Node] = pMakeSingleJAR3DModel(Search,Param,Prior,loopType);

  if isempty(Node),
%    mkdir([MotifLibraryPath filesep 'trouble']);
%    movefile([MotifLibraryPath filesep MotifName '.mat'],[MotifLibraryPath filesep 'trouble' filesep MotifName '.mat']); 
    fprintf('@@@@@@@@@@@@ pCalculateTransitions: Motif %s could not be modeled\n', MotifName);
  else
    [Text,T3,T4,T5] = xFASTACandidates(Search.File,Search,1,MotifName);

    % ---------- extract sequences counts for transition matricies
  
    [cWW, initial, transition, core, initialTransition] = pExtractSequenceCounts(Text);
    cWWCounts = cWWCounts + cWW;
    initialCounts = initialCounts + initial;
    coreCounts = coreCounts + core;
    transitionCounts = transitionCounts + transition;
    initialTransitionCounts = initialTransitionCounts + initialTransition;
  end
end

cWW_M = cWWCounts/sum(cWWCounts); 
initial_M = initialCounts/sum(initialCounts);
core_M = coreCounts/sum(coreCounts);
transition_M(1,:) = transitionCounts(1,:)/sum(transitionCounts(1,:));
transition_M(2,:) = transitionCounts(2,:)/sum(transitionCounts(2,:));
transition_M(3,:) = transitionCounts(3,:)/sum(transitionCounts(3,:));
transition_M(4,:) = transitionCounts(4,:)/sum(transitionCounts(4,:));
initialTransition_M(1,:) = initialTransitionCounts(1,:)/sum(initialTransitionCounts(1,:));
initialTransition_M(2,:) = initialTransitionCounts(2,:)/sum(initialTransitionCounts(2,:));
initialTransition_M(3,:) = initialTransitionCounts(3,:)/sum(initialTransitionCounts(3,:));
initialTransition_M(4,:) = initialTransitionCounts(4,:)/sum(initialTransitionCounts(4,:));

save(TransitionFile, 'cWW_M', 'initial_M', 'transition_M');

fprintf('pCalculateTransitions:  cWW basepair probabilities are\n');
cWW_M

fprintf('pCalculateTransitions:  initial distribution is\n');
initial_M

fprintf('pCalculateTransitions:  transition matrix is\n');
transition_M

