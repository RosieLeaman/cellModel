function [vertices,angles] = findVerticesNewMaterialCircle(centre,numVertices,radiusOrAmount,amountMaterial)

if radiusOrAmount == 0
    % we passed the amount, so we find the radius by solving for what
    % radius will give that area
    radius = sqrt(amountMaterial/pi);
else
    % we passed the radius directly (for testing purposes, so just use
    % that)
    radius = amountMaterial;
end

vertices = zeros(numVertices,2);
angles = zeros(numVertices,1);

index = 0;
while index < numVertices
    
    theta = 0 + (2*pi*index)/numVertices;
    angles(index+1) = theta;
    
    vertices(index+1,1) = centre(1)+radius*cos(theta);
    vertices(index+1,2) = centre(2)+radius*sin(theta);
    
    index = index + 1;   
end

