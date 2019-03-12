function [newPoints,errors] = surfaceTensionPolygon(points,dt)

% initialise newPoints
newPoints = zeros(size(points));

% we need to know how many vertices there are

vertices = size(points,1);

% we only need to run one iteration


% we should fit the spline here, as there's no point doing it every
% time

% but there is a special case when ii=1, so deal with that first

% matlab doesn't seem to like having the singularity at the end of
% our interval, which is the opposite to what I expected
% so for ii=1 rearrange our points and re-find the spline so that
% it's not at the end

halfway = floor(vertices/2);
points2halfway = [points(halfway:end,:);points(1:halfway,:)];

[normals2,~,splineX2,splineY2] = findTangentSpline(points2halfway,1);

points2 = [points;points(1,:)];

[normals,~,splineX,splineY] = findTangentSpline(points2,1);

errors = zeros(vertices,1);

% iterate over each point in the polygon
for ii = 1:size(points,1)

    % find the flow strength
    
    % (halfway/2)
    if ii < 5
        % because we shifted the points round to get the second spline
        % the index for the new normals in normals2 is not the same as
        % the normal index but is shifted
        correctNormalIndex = size(points2halfway,1)-halfway+ii;
        un = calculateIntegral(points(ii,:),normals2(correctNormalIndex,:),splineX2,splineY2);
    
    elseif size(points,1) - ii < 5
        % is this correctNormalIndex correct???? CHECK
        correctNormalIndex = ii - halfway + 1;
        un = calculateIntegral(points(ii,:),normals2(correctNormalIndex,:),splineX2,splineY2);            
    
    else        
        un = calculateIntegral(points(ii,:),normals(ii,:),splineX,splineY);
    end

    % this has to be applied in the direction of the outward pointing
    % normal

    flow = -un.*normals(ii,:);

    % move that point by how much

    newPoints(ii,:) = points(ii,:) + flow*dt;
    
    errors(ii) = findDist(points(ii,:),newPoints(ii,:));

end

figure;plot(errors)

end

function result = calculateIntegral(r0,r0normal,splineX,splineY)

% we will be returning un

% we want to use integral. So we need a function handle.

fun = @(t)calcIntegrandVectorised(t,r0,r0normal,splineX,splineY);

result = integral(fun,0,1);

end