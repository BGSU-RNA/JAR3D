% pAlignSequencesToMotifGroup mimics what JAR3D should soon be able to do,
% using the Java programs to align a set of sequences to an SCFG model and
% then using Python programs to show how those sequences correspond to the
% instances from 3D

if ~exist('loopType'),
  disp('Please specify a loop type, for example loopType = ''IL'';')
  return
end

OutputPath = ['Correspondences' filesep];

Filenames = dir(['MotifLibrary' filesep]);

keep = [];                               % of all models, which to keep

for m = 1:length(Filenames),
  if (length(Filenames(m).name) > 2),
    if (Filenames(m).name(1:2) == loopType),
      keep(m) = 1;
      Filenames(m).name = strrep(Filenames(m).name,'.mat','');
    end
  end 
end

Filenames = Filenames(find(keep));

for m = 1:length(Filenames),
  MotifName = Filenames(m).name;
  load(['MotifLibrary' filesep MotifName '.mat']);
  
  MN = Filenames(m).name;
  FN = ['MotifLibrary' filesep MN '.mat'];

  FastaFile = [pwd filesep OutputPath MotifName '.fasta'];
  ModelFile = [pwd filesep OutputPath MotifName '_model.txt'];

  [T,Node,Text] = pModelCorrespondences(MotifName,Search);

  fid = fopen(FastaFile,'w');
  for t = 1:length(Text),
    fprintf(fid,'%s\n',Text{t});
  end
  fclose(fid);

  Text = pNodeToSCFGModelText(Node,5);

  fid = fopen(ModelFile,'w');
  for i = 1:length(Text),
    fprintf(fid,'%s\n', Text{i});
  end
  fclose(fid);

  corresp = JAR3DMatlab.ModelCorrespondences(FastaFile,ModelFile,length(Text)+2);

  T = [T cell(corresp)];

  fid = fopen([OutputPath MotifName '_correspondences.txt'],'w');
  for r = 1:length(T),
    fprintf(fid,'%s\n',T{r});
  end
  fclose(fid);

%  disp('pAllModelCorrespondences: paused');

end

FN = 'IL_012';
ModelName = FN;
load(['MotifLibrary' filesep FN]);
FastaFile = 'C:/cygwin/home/zirbel/JAR3D/sequences/IL_018_13_cWW-tSH-tHH-cSH-tWH-tHS-cWW.fasta';
ModelFile = 'C:/cygwin/home/zirbel/JAR3D/models/IL_018_13_cWW-tSH-tHH-cSH-tWH-tHS-cWW.txt';
