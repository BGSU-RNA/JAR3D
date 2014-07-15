
function [void] = zStackedBar(bins,heights,color,orientation)

if nargin < 3,
  color = [];
end

if isempty(color),
  colormap('default');
  mymap = colormap;
  color = mymap([8 24 32 36 40 44 48:64],:);
end

if nargin < 4,
  orientation = 'h';
end

if lower(orientation) == 'h',

for i = 1:length(bins)-1,
  for j = 1:(length(heights(1,:))-1),
    B = heights(i,j);
    T = heights(i,j+1);
    L = bins(i);
    R = bins(i+1);
    jj = min(j,length(color));
    if color(jj) == 'o',
      patch([L R R L L], [B B T T B], [255  	165  	0]/255);
    elseif length(color(1,:)) == 3,
      patch([L R R L L], [B B T T B], color(jj,:));
    else
      patch([L R R L L], [B B T T B], color(jj));
    end
    hold on
  end
end

else

for i = 1:length(bins)-1,
  for j = 1:(length(heights(1,:))-1),
    B = heights(i,j);
    T = heights(i,j+1);
    L = bins(i);
    R = bins(i+1);
    jj = min(j,length(color));
    if color(jj) == 'o',
      patch([B B T T B], [L R R L L], [255  	165  	0]/255);
    elseif length(color(1,:)) == 3,
      patch([L R R L L], [B B T T B], color(jj,:));
    else
      patch([L R R L L], [B B T T B], color(jj));
    end
    hold on
  end
end

end
