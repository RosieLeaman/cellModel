function result = surfaceTensionTest2(numVertices,a,b)

result = 0;

% for testing
format long

% initialise an image
im = {};
idx = 1; %idx is the index of the current image in our future tiff file
xmax = 10;
ymax = 10;

format long
% if there are no vertices input we just use a default value
if nargin < 1
    numVertices = 1000;
end

if nargin <= 1
    % no ellipse parameters specified, use defaults
    a = 4; % x axis stretch
    b = 2; % y axis stretch
end

% make the points
[points,angles] = findVerticesNewMaterialEllipse([0,0],numVertices,a,b);
%[points,angles] = findVerticesNewMaterialEllipseWithError([0,0],numVertices,a,b,10e-4);

newPoints = points; % we initialise newPoints here as the current points

% plot the starting position of the points
im = addFrameToTiff(im,points,xmax,ymax,0,1);

% we will run the algorithm for so many time steps and show the output
time = 0;
maxTime = 0.02; % 0.3
dt = 0.01;

% here we store anything we would like to plot
uns = zeros(numVertices,1);
flows = zeros(numVertices,2);

disp('starting loop')
count = 0;
%while time < maxTime
    disp(['time is ',num2str(time),' out of ',num2str(maxTime)])

while count < 1
    count = count + 1;
    
    errors = zeros(1,size(points,1));
    
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

%     assert(mean(abs(-a*sin(x)./sqrt(a*a*sin(x).^2+b*b*cos(x).^2)-tangents(:,1)')) < 10e-14,'Estimated tangent x component not close to correct')
%     assert(mean(abs(b*cos(x)./sqrt(a*a*sin(x).^2+b*b*cos(x).^2)-tangents(:,2)')) < 10e-14,'Estimated tangent y component not close to correct')
    
    currentArea = polyarea(points(:,1),points(:,2));
    % we then iterate over each point in the circle
    for ii = 1:size(points,1)
    %for ii = 1

        % find the flow strength
        % we do this by calculating the surface tension integral
        % which gives the force in the normal direction
        
        % we must first calculate the rotation matrix which brings the
        % point ii to have normal (1,0) and tangent (0,1)
        
        % the correctedIndex is the index
        rotMatrix = [[tangents(ii,2),-tangents(ii,1)];[-normals(ii,2),normals(ii,1)]];
        invRotMatrix = [[normals(ii,1),tangents(ii,1)];[normals(ii,2),tangents(ii,2)]];
%         
        if count > 1000000 && ii <= 2
            un = calculateIntegral(points(ii,:),rotMatrix,splineX,splineY,1,0,a,b);
        else
            un = calculateIntegral(points(ii,:),rotMatrix,splineX,splineY,0,0,a,b);
        end
        %un = calculateIntegral(points(ii,:),rotMatrix,splineX,splineY,0,0,a,b);

%         if ii==1
%             disp('doing some checks')
%             t = linspace(0,1,10000);
%             
%             integrand = calcIntegrandVectorised(t,points(ii,:),rotMatrix,splineX,splineY,1,1,a,b);
%         end
        
        uns(ii) = un;
        
        % this has to be applied in the direction of the outward pointing
        % normal      
%         flow = -un.*normals(ii,:);
%         flows(ii,:) = flow;

        % we have the flow along the normal direction.
        
        newPoints(ii,:) = points(ii,:) + dt*un*normals(ii,:);
        flows(ii,:) =  dt*un*normals(ii,:);
        
        errors(ii) = findDist(points(ii,:),newPoints(ii,:));
        errors(ii) = norm(flows(ii,:));
      
    end

    % the error is just un really
    errors;
    calculatedError = mean(errors);
    %a = [min(errors),mean(errors),max(errors)];
    
    if count > 10000000 || time == maxTime
        figure; hold on;
        plot(points(:,1),points(:,2),'x-')
        %plot(newPoints(:,1),newPoints(:,2)','o-')
        xlim([-10,10]);ylim([-10,10]);
        title(count)
    end

    % move points to the new points
    points = newPoints;
    
    if mod(count,5) == 0
    %if 1
        im = addFrameToTiff(im,points,xmax,ymax,time,1);
    end
     
    time = time + dt;
    
    if 1
    %if count > 10000000
        figure;
        plot(angles,uns)
        xticks(0:pi/2:2*pi)
        xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
        ylabel('un');xlabel('angle')
        title(count)
        
        result = uns;

        figure; hold on;
        plot(angles,flows(:,1))
        plot(angles,flows(:,2))
        xticks(0:pi/2:2*pi)
        xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
        legend('x','y')
        ylabel('flows after rotation');xlabel('angle')
        grid minor;
        title(count)

    end
    
    % the best we ever do on the distance moved is 10e-10. So the absolute
    % best we can do for the squared distance change is 10e-5. This is also
    % an absolute best, it's possible we may do worse. Maybe 10e-4 would be
    % a good test then?
    % Note, we are currently doing worse than this. (About 10e-3 on
    % average)
    
    newArea = polyarea(points(:,1),points(:,2));
    count
    result = abs(1-(newArea/currentArea));
    
    result = uns;
    
    assert(abs(newArea - currentArea) < 10e-5,['Polygon changed area when calculating surface tension, was ',num2str(currentArea),' now ',num2str(newArea)]);
    
    
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

function result = calculateIntegral(r0,r0rotation,splineX,splineY,plotYes,normalYes,a,b)
% we will be returning un
% we want to use integral. So we need a function handle.

fun = @(t)calcIntegrandVectorised(t,r0,r0rotation,splineX,splineY,plotYes,normalYes,a,b);

result = integral(fun,0,1);

fun5 = @(z) (2*pi*a*(a^2*sin(2*pi*z).^2 + b^2*cos(2*pi*z).^2).^(1/2).*(cos(2*pi*z) - 1).*((a^2.*sin(2*pi*z).^2.*(b^2*sin(2*pi*z).^2 - a^2.*(cos(2*pi*z) - 1).^2))./(a^2.*sin(2*pi*z).^2 + b^2.*cos(2*pi*z).^2) - (b^2*cos(2*pi*z).^2.*(b^2.*sin(2*pi*z).^2 - a^2.*(cos(2*pi*z) - 1).^2))./(a^2.*sin(2*pi*z).^2 + b^2.*cos(2*pi*z).^2) + (4*a^2*b^2.*cos(2*pi*z).*sin(2*pi*z).^2.*(cos(2*pi*z) - 1))./(a^2.*sin(2*pi*z).^2 + b^2.*cos(2*pi*z).^2)))./(b^2.*sin(2*pi*z).^2 + a^2.*(cos(2*pi*z) - 1).^2).^2;


%actualIntegralX = integral(fun5,0,1)

end

function [image] = addFrameToTiff(image,points,xmax,ymax,time,closeFrameYes)
    % find the next index
    idx = numel(image)+1;
    
    % plot the image
    fig = figure;
    plot(points(:,1),points(:,2),'x');
    xlim([-xmax,xmax]);ylim([-ymax,ymax]);
    %xlim([3.6,4]);ylim([-1,1]);
    
    text(-xmax+1,-ymax+1,num2str(time))

    frame = getframe(fig);
    image{idx} = frame2im(frame);
    
    if closeFrameYes == 1
        close(fig)
    end
    
end