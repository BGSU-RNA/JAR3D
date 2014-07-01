% pJAR3DMaster runs JAR3D to make SCFG models, produce diagnostics

Release = ['HL' filesep '1.14'];
MotifLibraryLocation = 'C:\Users\zirbel\Documents\JAR3DMotifs\';
MotifLibraryPath = [MotifLibraryLocation filesep Release filesep 'motifs'];

OutputBase = pwd;                                                      % start Matlab in the directory containing IL and HL
JAR3Dpath = 'C:\Users\zirbel\Documents\GitHub\JAR3D\target\classes';   % location of class files for Java programs; needed for scoring sequences
javaaddpath(JAR3Dpath)

% -------------------------------------------------------------------- % run JAR3D on the selected motif group

if isempty(MotifLibraryPath),
	fprintf('No motif files found for %s\n',Release);
else
	pMakeSCFGModels(MotifLibraryPath,OutputBase,Release,1)

	pJAR3DDiagnostics(OutputBase,Release,Release,1)

	pGenerateRandomMotifSequences(OutputBase,Release)

	pJAR3DFalsePositiveStudy(OutputBase,Release,1)                         % 1 means to read FASTA files and parse and calculate edit distance; slow
	pJAR3DFalsePositiveStudy(OutputBase,Release,2)                         % 2 means to accumulate false positive data, which takes lots of RAM

	pSetModelSpecificCutoffs(OutputBase,Release,0)
	pJAR3DFalsePositiveStudy(OutputBase,Release,3)                         % 3 means to use model-specific cutoffs
end

