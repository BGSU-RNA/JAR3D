% pGetModelNames reads all _model.txt or .mat files in the specified folder

function [NumModels,MotifNames,ModNames,ModNums] = pGetModelNames(ModelPath,loopType)

% --------------------------------- read model filenames, select right ones

Filenames = dir(ModelPath);       % read motif library .mat files

keep = [];                               % of all models, which to keep

for m = 1:length(Filenames),
  N = Filenames(m).name;
  if (length(Filenames(m).name) > 2),
    if strcmp(N(1:2),loopType) && (~isempty(strfind(N,'_model.txt'))) || (~isempty(strfind(N,'.mat')) && (isempty(strfind(N,'lengthdist')))),
      keep(m) = 1;
    end
  end 
end

Filenames = Filenames(find(keep));

% --------------------------------- Make model file names

NumModels = length(Filenames);

if NumModels > 0,
	for m = 1:NumModels,
	  MotifNames{m} = strrep(Filenames(m).name,'_model.txt','');
	  ModNames{m} = Filenames(m).name;
	  ModNums{m} = ModNames{m}(4:8);        % extract 5-digit number
	end
else
	fprintf('No model files found in %s\n',ModelPath);
	MotifNames = {};
	ModNames = {};
	ModNums = {};
end
