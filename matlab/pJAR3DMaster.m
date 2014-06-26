% pJAR3DMaster runs JAR3D to make SCFG models, produce diagnostics

Release = ['IL' filesep '1.15'];
MotifLibraryLocation = 'C:\Users\zirbel\Documents\Motifs\';
MotifLibraryPath = '';

switch Release
case ['IL' filesep '1.4']
	MotifLibraryPath = [MotifLibraryLocation filesep 'IL_20130803_0657\mat'];
case ['HL' filesep '1.4']
	MotifLibraryPath = [MotifLibraryLocation filesep 'HL_20130803_1418\mat'];
case ['IL' filesep '1.5']
	MotifLibraryPath = [MotifLibraryLocation filesep 'IL_20130803_0657\mat'];
case ['HL' filesep '1.5']
	MotifLibraryPath = [MotifLibraryLocation filesep 'HL_20130803_1418\mat'];
case ['IL' filesep '1.6']
	MotifLibraryPath = [MotifLibraryLocation filesep 'IL_20130831_0706\mat'];
case ['HL' filesep '1.6']
	MotifLibraryPath = [MotifLibraryLocation filesep 'HL_20130831_1253\mat'];
case ['IL' filesep '1.7']
	MotifLibraryPath = [MotifLibraryLocation filesep 'IL_20130928_0607\mat'];
case ['HL' filesep '1.7']
	MotifLibraryPath = [MotifLibraryLocation filesep 'HL_20130928_1039\mat'];
case ['IL' filesep '1.8']
	MotifLibraryPath = [MotifLibraryLocation filesep 'IL_20131105_1114\mat'];
case ['HL' filesep '1.8']
	MotifLibraryPath = [MotifLibraryLocation filesep 'HL_20131103_1536\mat'];
case ['IL' filesep '1.9']
	MotifLibraryPath = [MotifLibraryLocation filesep '\mat'];
case ['HL' filesep '1.9']
	MotifLibraryPath = [MotifLibraryLocation filesep '\mat'];
case ['IL' filesep '1.10']
	MotifLibraryPath = [MotifLibraryLocation filesep '\mat'];
case ['HL' filesep '1.10']
	MotifLibraryPath = [MotifLibraryLocation filesep '\mat'];
case ['IL' filesep '1.11']
	MotifLibraryPath = [MotifLibraryLocation filesep '\mat'];
case ['HL' filesep '1.11']
	MotifLibraryPath = [MotifLibraryLocation filesep '\mat'];
case ['IL' filesep '1.12']
	MotifLibraryPath = [MotifLibraryLocation filesep '\mat'];
case ['HL' filesep '1.12']
	MotifLibraryPath = [MotifLibraryLocation filesep '\mat'];
case ['HL' filesep '1.13']
	MotifLibraryPath = [MotifLibraryLocation filesep 'HL_20140329_1054\mat'];
case ['IL' filesep '1.13']
	MotifLibraryPath = [MotifLibraryLocation filesep 'IL_20140329_0627\mat'];
case ['IL' filesep '1.14']
	MotifLibraryPath = [MotifLibraryLocation filesep '\mat'];
case ['HL' filesep '1.14']
	MotifLibraryPath = [MotifLibraryLocation filesep '\mat'];
case ['IL' filesep '1.15']
	MotifLibraryPath = [MotifLibraryLocation filesep 'IL_20140619_1135\mat'];
case ['HL' filesep '1.15']
	MotifLibraryPath = [MotifLibraryLocation filesep 'HL_20140614_2131\mat'];
end

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

