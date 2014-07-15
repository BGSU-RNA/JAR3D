% pJAR3DMaster runs JAR3D to make SCFG models, produce diagnostics

% -------------------------------------------------------------------- % choices for the user

Release = ['IL' filesep '1.3'];
Release = ['IL' filesep '1.12'];
Release = ['IL' filesep '1.14'];
Release = ['IL' filesep '1.13'];
Release = ['IL' filesep '1.2'];

MotifLibraryLocation = 'C:\Users\zirbel\Documents\JAR3DMotifs\';
OutputBase = pwd;                                                      % start Matlab in the directory containing IL and HL
JAR3Dpath = 'C:\Users\zirbel\Documents\GitHub\JAR3D\target\classes';   % location of class files for Java programs; needed for scoring sequences
Pythonpath = 'C:\Users\zirbel\Documents\GitHub\JAR3D\python';          % location of python programs in the JAR3D release

% -------------------------------------------------------------------- % run JAR3D on the selected motif group

MotifLibraryPath = [MotifLibraryLocation filesep Release filesep 'motifs'];
javaaddpath(JAR3Dpath)

pMakeSCFGModels(MotifLibraryPath,OutputBase,Release,1)
try
	system(['python ' Pythonpath filesep 'GroupToModelDiagnostic.py ' OutputBase filesep Release]);
catch
	fprintf('Not able to run the diagnostic in which sequences from 3D structures are aligned to their JAR3D models\n');
end
pJAR3DDiagnostics(OutputBase,Release,Release,1)

pGenerateRandomMotifSequences(OutputBase,Release)

pJAR3DFalsePositiveStudy(OutputBase,Release,1)                         % 1 means to read FASTA files and parse and calculate edit distance; slow
pJAR3DFalsePositiveStudy(OutputBase,Release,2)                         % 2 means to accumulate false positive data, which takes lots of RAM

pSetModelSpecificCutoffs(OutputBase,Release,0)
pJAR3DFalsePositiveStudy(OutputBase,Release,3)                         % 3 means to use model-specific cutoffs

