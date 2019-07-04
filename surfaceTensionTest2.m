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
    a = 2; % x axis stretch
    b = 2; % y axis stretch
end

[points,angles] = findVerticesNewMaterialEllipse([0,0],numVertices,a,b);
%[points,angles] = findVerticesNewMaterialEllipseWithError([0,0],numVertices,a,b,10e-4);

newPoints = points;

% plot the starting position of the points
im = addFrameToTiff(im,points,xmax,ymax,0,1);

% we will run the algorithm for so many time steps and show the output
time = 0;
maxTime = 0.02; % 0.3
dt = 0.01;

% store the error here, this is only calculated properly if we're only
% doing one time step otherwise it's not very meaningful
calculatedError = 0;

% here we store the flows at each vertex and the amount moved in the normal
% direction
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
    % we have to add the beginning of points to the end to give the
    % impression of a complete polygon as we don't store the repeated point
    
    %points2 = [points;points(1,:)];
    % alternative, wrap around our polygon three times
    points2 = [points;points;points];
    
    [splineX,splineY] = fitSpline(points2,parametrisation);
    [splineX2,splineY2] = fitSpline(points2,1);
    [wholeTangents,wholeNormals] = findTangentFromSplines(splineX.breaks,splineX,splineY,1);
    
    % the above gives us normals and tangents for each point going around
    % the polygon three times. Just select out the middle bit that we
    % actually want here
    normals = wholeNormals(size(points,1)+1:size(points,1)*2,:);
    tangents = wholeTangents(size(points,1)+1:size(points,1)*2,:);
    
%     figure; hold on;
%     for i=1:25:size(points,1)
%         plot(points(i,1),points(i,2),'bx')
%         
%         plot([points(i,1),points(i,1)+normals(i,1)],[points(i,2),points(i,2)+normals(i,2)],'r-o')
%         plot([points(i,1),points(i,1)+tangents(i,1)],[points(i,2),points(i,2)+tangents(i,2)],'k-o')
%     end
%     xlim([-10,10]);ylim([-10,10])
    
    xPoints = ppval(splineX,splineX.breaks(1:numVertices))';
    yPoints = ppval(splineY,splineY.breaks(1:numVertices))';

    % test the interpolation worked correctly
%     assert(mean(abs(xPoints-points(:,1))) < 10e-14,'Interpolated x points not close to actual points')
%     assert(mean(abs(yPoints-points(:,2))) < 10e-14,'Interpolated y points not close to actual points')
% %     
    % test the tangent is close to accurate
    % the analytical unit tangent should be for an ellipse
    % t = (-(a/b)*y,(b/a)*x). With unit tangent being t./||t||
    x = linspace(0,2*pi*(1-1/numVertices),numVertices);

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
    
%     assert(abs(newArea - currentArea) < 10e-5,['Polygon changed area when calculating surface tension, was ',num2str(currentArea),' now ',num2str(newArea)]);
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