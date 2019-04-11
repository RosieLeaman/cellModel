function visualiseSimple(model)

hold on;

% plotting of BAMlocs and regions needs to be separated. Same for LPS.

if numel(model.BAMlocs) > 0
    plot(model.BAMlocs(:,1),model.BAMlocs(:,2),'bd','linewidth',2)
end

if numel(model.proteinVertices) > 0
    for poly = 1:numel(model.proteinVertices)
    plot(model.proteinVertices{poly}(:,1),model.proteinVertices{poly}(:,2),'bx')
    end
end

if numel(model.LptDlocs) > 0
    plot(model.LptDlocs(:,1),model.LptDlocs(:,2),'rd','linewidth',2)
end

if numel(model.lpsVertices) > 0
    for poly = 1:numel(model.lpsVertices)
        plot(model.lpsVertices{poly}(:,1),model.lpsVertices{poly}(:,2),'rx')
    end
end
