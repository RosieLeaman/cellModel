function visualiseSimple(BAMlocs,proteinVertices,LptDlocs,lpsVertices,rightEdge,leftEdge)

hold on;

% plotting of BAMlocs and regions needs to be separated. Same for LPS.

if numel(BAMlocs) > 0
    plot(BAMlocs(:,1),BAMlocs(:,2),'bd','linewidth',2)
end

if numel(proteinVertices) > 0
    for poly = 1:numel(proteinVertices)
        plot(proteinVertices{poly}(:,1),proteinVertices{poly}(:,2),'bx')
    end
end

if numel(LptDlocs) > 0
    plot(LptDlocs(:,1),LptDlocs(:,2),'rd','linewidth',2)
end

if numel(lpsVertices) > 0
    for poly = 1:numel(lpsVertices)
        plot(lpsVertices{poly}(:,1),lpsVertices{poly}(:,2),'rx')
    end
end

% plot the left and right edges of bacteria
plot(rightEdge(:,1),rightEdge(:,2),'kx-')
plot(leftEdge(:,1),leftEdge(:,2),'kx-')