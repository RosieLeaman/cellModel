function [vertices,angles] = findVerticesNewMaterialEllipse(centre,numVertices,a,b)

vertices = zeros(numVertices,2);
angles = zeros(numVertices,1);

index = 0;
while index < numVertices
    
    theta = 0 + (2*pi*index)/numVertices;
    angles(index+1) = theta;
    
    vertices(index+1,1) = centre(1)+a*cos(theta);
    vertices(index+1,2) = centre(2)+b*sin(theta);
    
    index = index + 1;   
end