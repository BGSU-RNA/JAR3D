% pJAR3DMaster runs JAR3D to make SCFG models, produce diagnostics

% -------------------------------------------------------------------- % choices for the user

Release = ['HL' filesep '1.17'];
Release = ['IL' filesep '3.50'];

Input = ['HL_3.50_2021-10-13_17'];
Release = ['HL' filesep '3.50'];

Input = ['IL_3.50_2021-10-13_14'];
Release = ['IL' filesep '3.50'];

MotifLibraryLocation = 'C:\Users\zirbel\Documents\JAR3D\Motifs\';            % where to find the motif atlas release
OutputBase = 'C:\Users\zirbel\Documents\JAR3D\';                       % where to write JAR3D models
MatlabCodepath = 'C:\Users\zirbel\Documents\GitHub\JAR3D\Matlab';      % location of Matlab code
JAR3Dpath = 'C:\Users\zirbel\Documents\GitHub\JAR3D\target\classes';   % location of class files for Java programs; needed for scoring sequences
Pythonpath = 'C:\Users\zirbel\Documents\GitHub\JAR3D\python';          % location of python programs in the JAR3D release

% -------------------------------------------------------------------- % run JAR3D on the selected motif group

addpath(MatlabCodepath)
javaaddpath(JAR3Dpath)

MotifLibraryPath = [MotifLibraryLocation Input filesep 'mat'];

diary([OutputBase Release filesep 'JAR3D log' strrep(char(datetime),':','-') '.txt'])

pMakeSCFGModels(MotifLibraryPath,OutputBase,Release,1)

pGenerateRandomMotifSequences(OutputBase,Release,1)
pGenerateRandomMotifSequences(OutputBase,Release,2)
pGenerateRandomMotifSequences(OutputBase,Release,3)

pJAR3DFalsePositiveStudy(OutputBase,Release,1)  % 1 means to read FASTA files and parse and calculate edit distance; slow

pJAR3DFalsePositiveStudy(OutputBase,Release,2)  % 2 means to accumulate false positive data, which takes some RAM
pSetModelSpecificCutoffs(OutputBase,Release,0)  % Use results from FalsePositiveStudy Mode 2

pJAR3DFalsePositiveStudy(OutputBase,Release,3)  % 3 means to use model-specific cutoffs

pJAR3DDiagnostics(OutputBase,Release,Release,1)

try
	system(['python ' Pythonpath filesep 'GroupToModelDiagnostic.py ' OutputBase filesep Release]);
catch
	fprintf('Not able to run the diagnostic in which sequences from 3D structures are aligned to their JAR3D models\n');
end

