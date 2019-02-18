function visualiseSimple(model)

hold on;

if numel(model.BAMlocs) > 0
    for poly = 1:size(model.proteinVertices,3)
        plot(model.proteinVertices(:,1,poly),model.proteinVertices(:,2,poly),'bx')
    end

    plot(model.BAMlocs(:,1),model.BAMlocs(:,2),'bd','linewidth',2)
end

if numel(model.LptDlocs) > 0
    plot(model.LptDlocs(:,1),model.LptDlocs(:,2),'rd','linewidth',2)

    for poly = 1:size(model.lpsVertices,3)
        plot(model.lpsVertices(:,1,poly),model.lpsVertices(:,2,poly),'rx')
    end
end
