function calculatedError = surfaceTensionTest2(numVertices)

% initialise an image
im = {};
idx = 1; %idx is the index of the current image in our future tiff file
xmax = 10;
ymax = 10;

format long
% if there are no vertices input we just use a default value
if nargin < 1
    numVertices = 51;
end

a = 4; % x axis stretch
b = 2; % y axis stretch

[points,angles] = findVerticesNewMaterialEllipse([0,0],numVertices,a,b);

newPoints = points;

% plot the starting position of the points
im = addFrameToTiff(im,points,xmax,ymax,1);

% we will run the algorithm for so many time steps and show the output
time = 0;
maxTime = 0.6; % 0.3
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
    
    [splineX,splineY] = fitSpline(points2);
    [wholeTangents,wholeNormals] = findTangentFromSplines(splineX.breaks,splineX,splineY,1);
    
    % the above gives us normals and tangents for each point going around
    % the polygon three times. Just select out the middle bit that we
    % actually want here
    normals = wholeNormals(size(points,1)+1:size(points,1)*2,:);
    tangents = wholeTangents(size(points,1)+1:size(points,1)*2,:);
    
    xPoints = ppval(splineX,linspace(0,1-1/numVertices,numVertices))';
    yPoints = ppval(splineY,linspace(0,1-1/numVertices,numVertices))';

    % test the interpolation worked correctly
    assert(mean(abs(xPoints-points(:,1))) < 10e-14,'Interpolated x points not close to actual points')
    assert(mean(abs(yPoints-points(:,2))) < 10e-14,'Interpolated y points not close to actual points')
    
    % test the tangent is close to accurate
    % the analytical unit tangent should be for an ellipse
    % t = (-(a/b)*y,(b/a)*x). With unit tangent being t./||t||

    x = linspace(0,2*pi*(1-1/numVertices),numVertices);
    assert(mean(abs(-a*sin(x)./sqrt(a*a*sin(x).^2+b*b*cos(x).^2)-tangents(:,1)')) < 10e-14,'Estimated tangent x component not close to correct')
    assert(mean(abs(b*cos(x)./sqrt(a*a*sin(x).^2+b*b*cos(x).^2)-tangents(:,2)')) < 10e-14,'Estimated tangent y component not close to correct')

    % we then iterate over each point in the circle
    %for ii = 1:size(points,1)
    for ii = 1

        % find the flow strength
        % we do this by calculating the surface tension integral
        % which gives the force in the normal direction

        un = calculateIntegral(points(ii,:),normals(ii,:),splineX,splineY,0,0,a,b);

        if ii==1
            disp('doing some checks')
            t = linspace(0,1,10000);
            
            integrand = calcIntegrandVectorised(t,points(ii,:),normals(ii,:),splineX,splineY,1,1,a,b);
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

    % the error is just un really
    errors;
    calculatedError = mean(errors);

    points = newPoints;  
    
    if mod(count,2) == 0
        im = addFrameToTiff(im,points,xmax,ymax,0);
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


actualIntegralX = integral(fun5,0,1)
actualIntegralY = integral(fun5Y,0,1)

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