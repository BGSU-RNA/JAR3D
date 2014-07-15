% pMotifCollectionCharacteristics(GroupData) displays plots and data telling about motifs

function [void] = pMotifCollectionCharacteristics(GroupData)

Data = [cat(1,GroupData.NumNT) cat(1,GroupData.NumBasepairs) cat(1,GroupData.NumStacks) cat(1,GroupData.NumBPh) cat(1,GroupData.NumBR) cat(1,GroupData.NumInstances)];
Data = full(Data);

Names{1} = 'NumNT';
Names{2} = 'NumBasepairs';
Names{3} = 'NumStacks';
Names{4} = 'NumBPh';
Names{5} = 'NumBR';
Names{6} = 'NumInstances';

figure(1)
clf
for v = 1:6,
  subplot(2,3,v)
  n = hist(Data(:,v),0:max(Data(:,v)));
  hist(Data(:,v),0:max(Data(:,v)));
  axis([-0.5 max(Data(:,v))+0.5 0 1.1*max(n)]);
  xlabel(Names{v});
  ylabel('Number of motif groups');
end
subplot(2,3,1)

figure(2)
clf
c = 0;
for v = 1:5,
  rv = 0.5*(rand(size(Data(:,v)))-0.5);
  for w = (v+1):6,
    c = c + 1;
    subplot(3,5,c)
    rw = 0.5*(rand(size(Data(:,v)))-0.5);
    
    plot(Data(:,v)+rv,Data(:,w)+rw,'.');
    xlabel(Names{v});
    ylabel(Names{w});
  end
end
subplot(3,5,1)
