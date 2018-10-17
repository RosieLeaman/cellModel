%

% INPUTS:
% points, an nx2 vector where the columns are the points (x,y)

function result = lineIntegral(points,values)

% first we find the distance between these points

% add the starting point to the end
%points(end+1,:) = points(1,:);

distances = zeros(1,size(points,1));

for i=2:size(points,1)
    distances(i) = distances(i-1)+findDist(points(i,:),points(i-1,:));
end


% then we use trapz to calculate the integral

result = trapz(distances,values);
