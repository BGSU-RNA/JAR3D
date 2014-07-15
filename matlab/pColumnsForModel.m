% pColumnsForModel lists out what column each part of a JAR3D model should go into when displayed

function [T] = pColumnsForModel(Node,ModelName)

T = pColumnsForModelRecursive(Node,1);          % start with Node 1

for i = 1:length(T),
  T{i} = [ModelName '_' T{i} ' appears_in_column ' num2str(i)];
end
