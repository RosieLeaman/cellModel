function vertices = findVerticesNewMaterialCircle(centre,numVertices,amountMaterial)

radius = sqrt(amountMaterial/pi);

vertices = zeros(numVertices,2);

index = 0;

while index < numVertices
    
    theta = 0 + (2*pi*index)/numVertices;
    
    vertices(index+1,1) = centre(1)+radius*cos(theta);
    vertices(index+1,2) = centre(2)+radius*sin(theta);
    
    index = index + 1;   
end

