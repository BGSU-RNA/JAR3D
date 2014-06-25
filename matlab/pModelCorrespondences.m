% pModelCorrespondences() loads a motif group, makes a JAR3D model, writes
% out a fasta file, and writes out all of the correspondences

function [T,Node,Text] = pModelCorrespondences(ModelName,Search,OutputPath)

FN = ModelName;

JAR3D_path

[filePaths,Node] = pWriteSingleJAR3DModel(Search,FN,1,10);

T2 = pAlignMotifGroupToModel(Search,Node,FN,FN); 

[Text,T3,T4,T5] = xFASTACandidates(Search.File,Search,1,FN);

T6 = pColumnsForModel(Node,FN);

T = [T2 T3 T4 T5 T6];

if 0 > 1,
  FN = 'IL_01239.1.mat';
  ModelName = FN;
  load(['MotifLibrary' filesep FN]);
  FastaFile = ['C:\Documents and Settings\zirbel\My Documents\My Dropbox\BGSURNA\Motifs\Correspondences\' FN '.fasta'];
  ModelFile = ['C:\Documents and Settings\zirbel\My Documents\My Dropbox\BGSURNA\Motifs\Correspondences\' FN '.txt'];
end
