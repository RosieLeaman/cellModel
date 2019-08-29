function newPoints = surfaceTensionPolygon(points,dt,surfaceTensionStrength)

% initialise newPoints
newPoints = zeros(size(points));

% we need to know how many vertices there are

numVertices = size(points,1);

% we only need to run one iteration

% to get the tangents and normals at each point, we fit a spline to the
% polygon which can be differentiated to give us these

% to get nice estimates of the tangent at the 'beginning'/'end' of the
% polygon list, we extend the points list by tripling it. Then the
% central 'points' list has good estimates of the tangent
points2 = [points;points;points];

% we fit a spline to these points
[splineX,splineY] = fitSpline(points2);

% from the spline we calculate a tangent and normals at the polygon
% points. We actually only need the normals which are used to move the
% points along the normal direction
[wholeTangents,wholeNormals] = findTangentFromSplines(splineX.breaks,splineX,splineY,1);

% the above gives us normals and tangents for each point going around
% the polygon three times. Just select out the middle bit that we
% actually want here
normals = wholeNormals(numVertices+1:numVertices*2,:);
tangents = wholeTangents(numVertices+1:numVertices*2,:);

%allTangents(count,:,:) = tangents;

%     assert(mean(abs(-a*sin(x)./sqrt(a*a*sin(x).^2+b*b*cos(x).^2)-tangents(:,1)')) < 10e-14,'Estimated tangent x component not close to correct')
%     assert(mean(abs(b*cos(x)./sqrt(a*a*sin(x).^2+b*b*cos(x).^2)-tangents(:,2)')) < 10e-14,'Estimated tangent y component not close to correct')

%currentArea = polyarea(points(:,1),points(:,2));

% we then iterate over each point in the circle
for ii = 1:numVertices
%for ii = 1

    % find the flow strength
    % we do this by calculating the surface tension integral
    % which gives the force in the normal direction

    % we must first calculate the rotation matrix which brings the
    % point ii to have normal (1,0) and tangent (0,1)

    rotMatrix = [[tangents(ii,2),-tangents(ii,1)];[-normals(ii,2),normals(ii,1)]];
    %invRotMatrix = [[normals(ii,1),tangents(ii,1)];[normals(ii,2),tangents(ii,2)]];

    % calculate the integral, which gives us the size of the flow in
    % the normal direction at the current point, points(ii,:)
    un = calculateIntegral(points(ii,:),rotMatrix,splineX,splineY,0);

    try
        assert(~isnan(un),'oops')
    catch
        fid = fopen('errorlog.txt','wt');
        fprintf(fid, 'Error: integral returned NaN');
        fclose(fid);
        if numel(im) > 0
            saveLocation = 'testEllipse.tif';
            imwrite(im{1},saveLocation)
            for i=2:numel(im)
                imwrite(im{i},saveLocation,'WriteMode','append')
            end
        end

        % figure;
        % plot(points(:,1),points(:,2),'x-')

        save('test.mat')
        assert(1==0)
    end

    % apply the flow along the outward pointing normal direction

    newPoints(ii,:) = points(ii,:) + surfaceTensionStrength*dt*un*normals(ii,:);
    
end

end

function result = calculateIntegral(r0,r0rotation,splineX,splineY,plotYes)
% we will be returning un
% we want to use integral. So we need a function handle.

fun = @(t)calcIntegrandVectorised(t,r0,r0rotation,splineX,splineY,plotYes);

result = integral(fun,0,1);

try
    assert(~isnan(result),'Integral returned NaN');
catch
    r0
    r0rotation
    fun = @(t)calcIntegrandVectorised(t,r0,r0rotation,splineX,splineY,1);

    result = integral(fun,0,1);
end

end