% plots the state of the cell with filled polygons. Requires the model
% (struct) as input

function visualise(model)

figure; hold on;
for poly = 1:size(model.proteinVertices,3)
    h = fill(model.proteinVertices(:,1,poly),model.proteinVertices(:,2,poly),'b');
    set(h,'facealpha',.5)
end

if numel(model.lpsVertices) > 0
    for poly = 1:size(model.lpsVertices,3)
        h = fill(model.lpsVertices(:,1,poly),model.lpsVertices(:,2,poly),'r');
        set(h,'facealpha',.5)
    end
    plot(model.LptDlocs(:,1),model.LptDlocs(:,2),'d','linewidth',2)
end

plot(model.BAMlocs(:,1),model.BAMlocs(:,2),'x','linewidth',2)

   

