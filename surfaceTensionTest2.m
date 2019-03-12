function error = surfaceTensionTest2(vertices)

% get a circle
if nargin < 1
    vertices = 51;
end

points = zeros(vertices,2);
angles = zeros(vertices,1);

index = 0;
while index < vertices
    theta = 0 + (2*pi*index)/vertices;
    points(index+1,1) = cos(theta);
    points(index+1,2) = sin(theta);
    angles(index+1) = theta;
    index = index + 1;
end

% for spout
%points(18,:) = 1.5*points(18,:);

newPoints = points;

% for a certain length of time

time = 0;
maxTime = 0.3;
dt = 0.01;

% errors!

error = 0;

figure;
plot(points(:,1),points(:,2),'x-');
xlim([-2.2,2.2]);ylim([-2.2,2.2]);
hold on;

%pause();

count = 0;
idx = 1;

disp('starting loop')

uns = zeros(vertices,1);

%while time < maxTime
while count < 1
    count = count + 1;
    
    errors = zeros(1,size(points,1));
    
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
    
    % iterate over each point in the circle
    for ii = 1:size(points,1)

        % find the flow strength
        
        if ii < (halfway/2)
            % because we shifted the points round to get the second spline
            % the index for the new normals in normals2 is not the same as
            % the normal index but is shifted
            correctNormalIndex = size(points2halfway,1)-halfway+ii;
            un = calculateIntegral(points(ii,:),normals2(correctNormalIndex,:),splineX2,splineY2);
        elseif size(points,1) - ii < (halfway/2)
            % is this correctNormalIndex correct???? CHECK
            correctNormalIndex = ii - halfway + 1;
            un = calculateIntegral(points(ii,:),normals2(correctNormalIndex,:),splineX2,splineY2);            
        else        
            un = calculateIntegral(points(ii,:),normals(ii,:),splineX,splineY);
            if ii==2
                t = linspace(0,1,500);
                integrand = calcIntegrandVectorised(t,points(ii,:),normals(ii,:),splineX,splineY);
                figure;
                plot(t,integrand);
                title('integrand point2')
            end
        end
        
        uns(ii) = un;
        
        % this has to be applied in the direction of the outward pointing
        % normal
        
        flow = -un.*normals(ii,:);

        % move that point by how much
        
        newPoints(ii,:) = points(ii,:) + flow*dt;
        
        errors(ii) = findDist(points(ii,:),newPoints(ii,:));
      
    end
    
    % the error is just un actually
    error = mean(errors);
    
    points = newPoints;  
    
    
    if mod(count,3) == 0
        fig = figure;
        plot(points(:,1),points(:,2),'x');
        xlim([-2.2,2.2]);ylim([-2.2,2.2]);
        
        frame = getframe(fig);
        im{idx} = frame2im(frame);
        idx = idx + 1;
    end
    
     
    time = time + dt;
    
%     figure;
%     plot(angles,uns)
%     
    figure;
    plot(angles,errors,'x-')
    
%     figure;
%     plot(angles,abs(cos(angles)-normals(1:end-1,1)))
%     title('x component')
%     
%     figure;
%     plot(angles,abs(sin(angles)-normals(1:end-1,2)))
%     title('y component')

end

% saveLocation = 'testCircle.tif';
% imwrite(im{1},saveLocation)
% for i=2:numel(im)
%     imwrite(im{i},saveLocation,'WriteMode','append')
% end

end

function result = calculateIntegral(r0,r0normal,splineX,splineY)

% we will be returning un

% we want to use integral. So we need a function handle.

fun = @(t)calcIntegrandVectorised(t,r0,r0normal,splineX,splineY);

result = integral(fun,0,1);

end