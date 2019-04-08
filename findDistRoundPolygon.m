function dist = findDistRoundPolygon(xIndex,yIndex,polygon)

x = polygon(xIndex,:);
y = polygon(yIndex,:);

% is it faster to go left or right around the polygon?

halfway = size(polygon,1)/2;

% 1 2    3456 7    89
% -xIndex----yIndex-- 9 total halfway = 4.5.
% yIndex - xIndex = 5 > halfway
% faster to go from yIndex to xIndex than xIndex to yIndex

if yIndex > xIndex
    if yIndex - xIndex <= halfway
        startIndex = xIndex;
        endIndex = yIndex;
    else
        startIndex = yIndex;
        endIndex = xIndex;
    end
else
    if xIndex - yIndex <= halfway
        startIndex = yIndex;
        endIndex = xIndex;
    else
        startIndex = xIndex;
        endIndex = yIndex;
    end  
end

% now we go from the startIndex, to endIndex now we know the order

dist = 0;

while startIndex ~= endIndex
    % add the distance to the next index
    if startIndex < size(polygon,1)
        dist = dist + findDist(polygon(startIndex,:),polygon(startIndex+1,:));
    else
        dist = dist + findDist(polygon(startIndex,:),polygon(1,:));
    end
    
    % shuffle startIndex up one
    startIndex = startIndex + 1;
    if startIndex > size(polygon,1)
        startIndex = 1;
    end
end