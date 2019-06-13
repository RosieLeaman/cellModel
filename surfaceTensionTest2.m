function error = surfaceTensionTest2(numVertices)

% initialise an image
im = {};
idx = 1; %idx is the index of the current image in our future tiff file
xmax = 5;
ymax = 5;

format long
% if there are no vertices input we just use a default value
if nargin < 1
    numVertices = 51;
end

a = 2; % x axis stretch
b = 2; % y axis stretch

[points,angles] = findVerticesNewMaterialEllipse([0,0],numVertices,a,b);

newPoints = points;

% plot the starting position of the points
im = addFrameToTiff(im,points,xmax,ymax,1);

% we will run the algorithm for so many time steps and show the output
time = 0;
maxTime = 0.1; % 0.3
dt = 0.01;

% store the error here, this is only calculated properly if we're only
% doing one time step otherwise it's not very meaningful
error = 0;

% here we store the flows at each vertex and the amount moved in the normal
% direction
uns = zeros(numVertices,1);
flows = zeros(numVertices,2);

disp('starting loop')
%while time < maxTime
count = 0;
while count < 1
    count = count + 1;
    
    errors = zeros(1,size(points,1));
    
    % to get the tangents and normals at each point, we fit a spline to the
    % polygon which can be differentiated to give us these
    % we have to add the beginning of points to the end to give the
    % impression of a complete polygon as we don't store the repeated point
    
    %points2 = [points;points(1,:)];
    % alternative, wrap around our polygon three times
    points2 = [points;points;points];
    
    [normals,tangents,splineX,splineY] = findTangentSpline(points2,1);
    
    xPoints = ppval(splineX,linspace(0,1-1/numVertices,numVertices))';
    yPoints = ppval(splineY,linspace(0,1-1/numVertices,numVertices))';

    % test the interpolation worked correctly
    assert(mean(abs(xPoints-points(:,1))) < 10e-14,'Interpolated x points not close to actual points')
    assert(mean(abs(yPoints-points(:,2))) < 10e-14,'Interpolated y points not close to actual points')
    
    % test the tangent is close to accurate
    % the analytical unit tangent should be for an ellipse
    % (-(a/b)*y,(b/a)*x)/||t||
    
    correctTangent = zeros(size(points));
    for i=1:size(points,1)
        correctTangent(i,1) = -(a/b)*points(i,2);
        correctTangent(i,2) = (b/a)*points(i,1);
        tangentNorm = norm(correctTangent(i,:));
        
        correctTangent(i,:) = correctTangent(i,:)./tangentNorm;
    end
    
    x = linspace(0,2*pi*(1-1/numVertices),1000);
    assert(mean(abs(-sin(x)-correctTangent(:,1)')) < 10e-14,'Estimated tangent x component not close to correct')
    assert(mean(abs(cos(x)-correctTangent(:,2)')) < 10e-14,'Estimated tangent y component not close to correct')

    break
    
%     angles2 = [angles;2*pi];
%     angles2 = angles2./(2*pi);
%     splineX = spline(angles2,points2(:,1));
%     splineY = spline(angles2,points2(:,2));
    
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
    
    halfway = floor(numVertices/2);
    points2halfway = [points(halfway:end,:);points(1:halfway,:)];
%     [normals2,~,splineX2,splineY2] = findTangentSpline(points2halfway,1);
%     

    %figure;plot(ppx,ppy,'kx-')

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
    %for ii = 1:size(points,1)
    for ii = 1

        % find the flow strength
        % we do this by calculating the surface tension integral
        % which gives the force in the normal direction
        
        % first we deal with points where we use the second spline
        if ii < (halfway/2)
            % because we shifted the points round to get the second spline
            % the index for the new normals in normals2 is not the same as
            % the normal index but is shifted
            correctNormalIndex = size(points2halfway,1)-halfway+ii;
            %un = calculateIntegral(points(ii,:),normals2(correctNormalIndex,:),splineX2,splineY2,0,0,a,b);
        elseif size(points,1) - ii < (halfway/2)
            % again we have to shift the indices
            correctNormalIndex = ii - halfway + 1;
            un = calculateIntegral(points(ii,:),normals2(correctNormalIndex,:),splineX2,splineY2,0,0,a,b);            
        else   
            un = calculateIntegral(points(ii,:),normals(ii,:),splineX,splineY,0,0,a,b);
        end
        
        un = calculateIntegral(points(ii,:),normals(ii,:),splineX,splineY,0,1,a,b);
        %un = calculateIntegral(points(ii,:),[1,0],splineX,splineY,1,0);

        if ii==1
            disp('doing some checks')
            correctNormalIndex = size(points2halfway,1)-halfway+ii;
            t = linspace(0,1,10000);
            integrand = calcIntegrandVectorised(t,points(ii,:),normals(ii,:),splineX,splineY,1,1,a,b);
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

if numel(im) > 0
    saveLocation = 'testEllipse.tif';
    imwrite(im{1},saveLocation)
    for i=2:numel(im)
        imwrite(im{i},saveLocation,'WriteMode','append')
    end
end

end

function result = calculateIntegral(r0,r0normal,splineX,splineY,plotYes,normalYes,a,b)
% we will be returning un
% we want to use integral. So we need a function handle.

fun = @(t)calcIntegrandVectorised(t,r0,r0normal,splineX,splineY,plotYes,normalYes,a,b);

result = integral(fun,0,1)

fun5 = @(z) (2*pi*a*(a^2*sin(2*pi*z).^2 + b^2*cos(2*pi*z).^2).^(1/2).*(cos(2*pi*z) - 1).*((a^2.*sin(2*pi*z).^2.*(b^2*sin(2*pi*z).^2 - a^2.*(cos(2*pi*z) - 1).^2))./(a^2.*sin(2*pi*z).^2 + b^2.*cos(2*pi*z).^2) - (b^2*cos(2*pi*z).^2.*(b^2.*sin(2*pi*z).^2 - a^2.*(cos(2*pi*z) - 1).^2))./(a^2.*sin(2*pi*z).^2 + b^2.*cos(2*pi*z).^2) + (4*a^2*b^2.*cos(2*pi*z).*sin(2*pi*z).^2.*(cos(2*pi*z) - 1))./(a^2.*sin(2*pi*z).^2 + b^2.*cos(2*pi*z).^2)))./(b^2.*sin(2*pi*z).^2 + a^2.*(cos(2*pi*z) - 1).^2).^2;

actualIntegral = integral(fun5,0,1)

fun5Alt = @(t) (a*(a^2*sin(t).^2 + b^2*cos(t).^2).^(1/2).*(cos(t) - 1).*((a^2.*sin(t).^2.*(b^2*sin(t).^2 - a^2.*(cos(t) - 1).^2))./(a^2.*sin(t).^2 + b^2.*cos(t).^2) - (b^2*cos(t).^2.*(b^2.*sin(t).^2 - a^2.*(cos(t) - 1).^2))./(a^2.*sin(t).^2 + b^2.*cos(t).^2) + (4*a^2*b^2.*cos(t).*sin(t).^2.*(cos(t) - 1))./(a^2.*sin(t).^2 + b^2.*cos(t).^2)))./(b^2.*sin(t).^2 + a^2.*(cos(t) - 1).^2).^2;

actualIntegral = integral(fun5Alt,0,2*pi)

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

function [image] = addFrameToTiff(image,points,xmax,ymax,closeFrameYes)
    % find the next index
    idx = numel(image)+1;
    
    % plot the image
    fig = figure;
    plot(points(:,1),points(:,2),'x');
    xlim([-xmax,xmax]);ylim([-ymax,ymax]);

    frame = getframe(fig);
    image{idx} = frame2im(frame);
    
    if closeFrameYes == 1
        close(fig)
    end
    
end