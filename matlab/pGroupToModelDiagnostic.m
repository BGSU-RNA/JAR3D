% pGroupToModelDiagnostic runs a python program which shows how the sequences of 3D instances in each motif group correspond to the conserved positions in their motif group, and shows how JAR3D aligns the sequences to the motif group.

% function [void] = pGroupToModelDiagnostic(OutputBase,Release,PythonDirectory)

Release = ['IL' filesep '1.13'];
MotifLibraryLocation = 'C:\Users\zirbel\Documents\JAR3DMotifs\';
MotifLibraryPath = [MotifLibraryLocation filesep Release filesep 'motifs'];

OutputBase = pwd;                                                      % start Matlab in the directory containing IL and HL
JAR3Dpath = 'C:\Users\zirbel\Documents\GitHub\JAR3D\target\classes';   % location of class files for Java programs; needed for scoring sequences
javaaddpath(JAR3Dpath)
Pythonpath = 'C:\Users\zirbel\Documents\GitHub\JAR3D\python';

system(['python ' Pythonpath filesep 'GroupToModelDiagnostic.py ' OutputBase filesep Release])

