function error = surfaceTensionTest2(vertices)

format long
% if there are no vertices input we just use test points
if nargin < 1
    vertices = 51;
end

% these points are an ellipse
points = zeros(vertices,2);
angles = zeros(vertices,1);

a = 2;
b = 2;

index = 0;
while index < vertices
    theta = 0 + (2*pi*index)/vertices;
    points(index+1,1) = a*cos(theta);
    points(index+1,2) = b*sin(theta);
    angles(index+1) = theta;
    index = index + 1;
end

% for spout
%points(18,:) = 1.5*points(18,:);

newPoints = points;

% plot the starting position of the points
idx = 1; %idx is the index of the current image in our future tiff file

fig = figure;
plot(points(:,1),points(:,2),'x');
xlim([-20,20]);ylim([-20,20]);

frame = getframe(fig);
im{idx} = frame2im(frame);
idx = idx + 1;
%close(fig)

% we will run the algorithm for so many time steps and show the output
time = 0;
maxTime = 0.1; % 0.3
dt = 0.01;

% store the error here, this is only calculated properly if we're only
% doing one time step otherwise it's not very meaningful
error = 0;
count = 0;

disp('starting loop')

% here we store the flows at each vertex and the amount moved in the normal
% direction
uns = zeros(vertices,1);
flows = zeros(vertices,2);

%while time < maxTime
while count < 1
    count = count + 1;
    
    errors = zeros(1,size(points,1));
    
    % to get the tangents and normals at each point, we fit a spline to the
    % polygon which can be differentiated to give us these
    % we have to add the beginning of points to the end to give the
    % impression of a complete polygon as we don't store the repeated point
    
    points2 = [points;points(1,:)];
    [normals,~,splineX,splineY,totalDistAround] = findTangentSpline(points2,1);
    
    % matlab doesn't seem to like having the singularity at the end of
    % our interval, which is the opposite to what I expected
    % so for ii=1 rearrange our points and re-find the spline so that
    % it's not at the end
    
    % the spline fit is bad at the ends of the interval, so for points near
    % the beginning/end of the polygon vertex list we actually use an
    % alternative spline, spline2 which begins halfway around the polygon
    % but otherwise uses the same set of points. Spline2 estimates points
    % near the middle of the polygon vertex list poorly, but should be good
    % at the end points.
    
    halfway = floor(vertices/2);
    points2halfway = [points(halfway:end,:);points(1:halfway,:)];
    [normals2,~,splineX2,splineY2,totalDistAround2] = findTangentSpline(points2halfway,1);
    
    if abs(totalDistAround - totalDistAround2) > 10e-10
        warning('The distance around polygon was different in both cases.')
        disp(totalDistAround)
        disp(totalDistAround2)
    end
    
    % we plot the radius of curvature of the polygon
    % we have to re-do the spline as findTangentSpline 
    % does not return the second derivatives which are needed in the
    % formula
    
    % first get the distances apart from each other
    t = zeros(size(points2,1),1);
    for i=2:numel(t)
        t(i) = t(i-1) + findDist(points2(i,:),points2(i-1,:));
    end

    % normalise to [0,1]
    t = t./(max(t));
    
    ppvalXdiff = fnder(splineX);
    ppvalYdiff = fnder(splineY);
    ppvalXdiff2 = fnder(ppvalXdiff);
    ppvalYdiff2 = fnder(ppvalYdiff);

    xdiffs = ppval(ppvalXdiff,t);
    ydiffs = ppval(ppvalYdiff,t);
    xdiffs2 = ppval(ppvalXdiff2,t);
    ydiffs2 = ppval(ppvalYdiff2,t);
    
    R = zeros(size(points2,1),1);

    for i=1:numel(xdiffs)
        R(i) = (xdiffs(i)*ydiffs2(i)-ydiffs(i)*xdiffs2(i))/(((xdiffs(i))^2+(ydiffs(i))^2)^(3/2));
    end
    
%     figure;
%     plot(R,'linewidth',2)
%     axis tight
%     grid on
%     set(gca,'FontSize',16)
%     title('radius of curvature')
    
    % we then iterate over each point in the circle
    for ii = 1:size(points,1)
    %for ii = 1

        % find the flow strength
        % we do this by calculating the surface tension integral
        % which gives the force in the normal direction
        
        % first we deal with points where we use the second spline
        if ii < (halfway/2)
            % because we shifted the points round to get the second spline
            % the index for the new normals in normals2 is not the same as
            % the normal index but is shifted
            correctNormalIndex = size(points2halfway,1)-halfway+ii;
            un = calculateIntegral(points(ii,:),normals2(correctNormalIndex,:),splineX2,splineY2,totalDistAround,0,0);
        elseif size(points,1) - ii < (halfway/2)
            % again we have to shift the indices
            correctNormalIndex = ii - halfway + 1;
            un = calculateIntegral(points(ii,:),normals2(correctNormalIndex,:),splineX2,splineY2,totalDistAround,0,0);            
        else   
            un = calculateIntegral(points(ii,:),normals(ii,:),splineX,splineY,totalDistAround,0,0);
        end
        
        un = calculateIntegral(points(ii,:),normals(ii,:),splineX,splineY,totalDistAround,0,0);

        if ii==1
            disp('doing some checks')
            correctNormalIndex = size(points2halfway,1)-halfway+ii;
            t = linspace(0,2*pi,10000);
            integrand = calcIntegrandVectorised(t,points(ii,:),normals2(correctNormalIndex,:),splineX,splineY,1,0);
%             figure;
%             plot(t,integrand);
%             title(['integrand point (',num2str(points(ii,1)),',',num2str(points(ii,2)),')'])

        end
        
        uns(ii) = un;
        
        % this has to be applied in the direction of the outward pointing
        % normal      
        flow = -un.*normals(ii,:);
        flows(ii,:) = flow;

        % move that point by how much
        
        newPoints(ii,:) = points(ii,:) + flow*dt;
        
        errors(ii) = findDist(points(ii,:),newPoints(ii,:));
      
    end

    % the error is just un actually
    errors;
    error = mean(errors);
    
    
    points = newPoints;  
    
    if mod(count,10) == 0
    %if 1
        fig = figure;
        plot(points(:,1),points(:,2),'x');
        xlim([-20,20]);ylim([-20,20]);
        
        frame = getframe(fig);
        im{idx} = frame2im(frame);
        idx = idx + 1;
        %close(fig)
    end
    
     
    time = time + dt;
    
%     figure;
%     plot(angles,uns)
%     ylabel('un');xlabel('angle')
%     
%     figure;
%     plot(angles,errors,'x-')
%     
%     figure;
%     plot(angles,abs(cos(angles)-normals(1:end-1,1)))
%     title('x component')
%     
%     figure;
%     plot(angles,abs(sin(angles)-normals(1:end-1,2)))
%     title('y component')

end

saveLocation = 'testEllipse.tif';
imwrite(im{1},saveLocation)
for i=2:numel(im)
    imwrite(im{i},saveLocation,'WriteMode','append')
end

end

function result = calculateIntegral(r0,r0normal,splineX,splineY,totalDistAround,plotYes,normalYes)

% we will be returning un
% we want to use integral. So we need a function handle.

fun = @(t)calcIntegrandVectorised(t,r0,r0normal,splineX,splineY,plotYes,normalYes);

result = integral(fun,0,totalDistAround);

end

% this is a helpful function
function previous = prev(x,maxi)
    if x > 1
        previous = x - 1;
    else
        previous = maxi;
    end  
end

function following = next(x,maxi)
    if x < maxi
        following = x + 1;
    else
        following = 1;
    end  
end

% produces the indices of the next N numbers going in a circle
function following = nextN(x,maxi,N)
    following = zeros(1,N);
    following(1) = x;
    for i=2:N
        following(i) = next(following(i-1),maxi);
    end

end

% shuffle the numbers in the 1xm list along by N, so that result(1) =
% list(N)
function result = shuffle(list,N)
    result = nextN(N,numel(list),numel(list));
end