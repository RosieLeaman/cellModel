% this plots just the protein and lps regions in a fancy way
% it goes to one level of depth

function visualisePretty(proteinVertices,lpsVertices)

figure;
hold on;

% protein regions are always bottom, plot these first
for i=1:numel(proteinVertices)
    fill(proteinVertices{i}(:,1),proteinVertices{i}(:,2),'b')
end

% then plot lps on top
for i=1:numel(lpsVertices)
    fill(lpsVertices{i}(:,1),lpsVertices{i}(:,2),'r')
end

% if the simulation runs long enough then protein regions may be spawned
% inside lps regions, add these on here
for i=1:numel(proteinVertices)
    inside = checkRegionInOther(proteinVertices{i},lpsVertices);
    if inside > 0
        fill(proteinVertices{i}(:,1),proteinVertices{i}(:,2),'b')
    end
end

% if we want we can go deeper
% check for lps inside protein which is inside lps
% this is probably a bit too deep tbh
for i=1:numel(lpsVertices)
    [inside,j] = checkRegionInOther(lpsVertices{i},proteinVertices);
    if inside > 0
        inside2 = checkRegionInOther(proteinVertices{j},lpsVertices);
        if inside2 > 0
            fill(lpsVertices{i}(:,1),lpsVertices{i}(:,2),'r')
        end
    end
end

end

function [inside, j] = checkRegionInOther(polygon,otherVertices)
    inside = 0;
    for j=1:numel(otherVertices)
        polygonInside = inpolygon(polygon(:,1),polygon(:,2),otherVertices{j}(:,1),otherVertices{j}(:,2));
        inside = sum(polygonInside);
        if inside > 0
            break
        end
    end
end
