% function takes in a cell array of polygons and returns an array with same
% length with 0 if [xq,yq] is not inside polygonCells{i} and 1 if it is 

function inside = pointInsidePolygon(xq,yq,polygonCells)

inside = zeros(1,numel(polygonCells));

for i=1:numel(polygonCells)
    in = inpolygon(xq,yq,polygonCells{i}(:,1),polygonCells{i}(:,2));
    inside(i) = in;
end