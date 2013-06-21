function meanPlot( cellDiagrams, mean )
%   meanPlot creates a common plot of a collection of diagrams and
%       there mean.

%   @param cellDiagrams is the cell array of diagrams used for sampling.
%   @param mean is the mean of the diagrams collected in cellDiagrams.


%   Fill in empty diagrams to prevent crash.
for i=1:size(cellDiagrams, 2)
    if size(cellDiagrams{i},1)==0&&size(cellDiagrams{i},2)==0
        cellDiagrams{i}=zeros(0,2);
    end
end

%   Plot first diagram and grab greatest x or y value for axis scaling.
plot(cellDiagrams{1,1}(:,1),cellDiagrams{1,1}(:,2),'.');
t = max(max(cellDiagrams{1,1}));

%   Loop to plot rest of diagram collection and find correct scaling value.
for cell = 2:size(cellDiagrams, 2)
    hold on
    plot(cellDiagrams{1,cell}(:,1),cellDiagrams{1,cell}(:,2),'.')
    s = max(max(cellDiagrams{1,cell}));
    
    if t < s
        t = s;
    end
end

%   Plot mean in different style for comparison. Also, update scaling
%       value and then scale axes.
hold all
plot(mean(:,1),mean(:,2),'ro')
s = max(max(mean));

if t < s
    t = s;
end

axis([0,t*1.1,0,t*1.1]);

%   Plot diagonal.
hold on
refline(1,0)
hold off
end

