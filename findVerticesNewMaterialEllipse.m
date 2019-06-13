function [points,angles] = findVerticesNewMaterialEllipse(centre,numVertices,xAxis,yAxis)

points = zeros(numVertices,2);
angles = zeros(numVertices,1);

index = 0;
while index < numVertices
    theta = 0 + 2*pi*index/numVertices;
    points(index+1,1) = xAxis*cos(theta);
    points(index+1,2) = yAxis*sin(theta);
    angles(index+1) = theta;
    index = index + 1;
end