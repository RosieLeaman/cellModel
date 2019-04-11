function in = checkLineInPolygon(startPoint,endPoint,polygon)

xq = linspace(startPoint(1),endPoint(1),10);
yq = linspace(startPoint(2),endPoint(2),10);

in = inpolygon(xq,yq,polygon(:,1),polygon(:,2));