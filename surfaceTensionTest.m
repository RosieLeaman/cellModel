function error = surfaceTensionTest(vertices)

% get a circle
if nargin < 1
    vertices = 51;
end

% radius = 1;
% amount = (vertices*radius/2)*sin(2*pi/vertices);
% points = findVerticesNewMaterial([0,0],vertices,amount);

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


newPoints = points;

% for a certain length of time

time = 0;
maxTime = 0.1;
dt = 0.01;

% errors!

error = 0;

figure;
plot(points(:,1),points(:,2),'x-');
xlim([-1.4,1.4]);ylim([-1.4,1.4]);
hold on;

%pause();

count = 0;
idx = 1;

disp('starting loop')

uns = zeros(vertices,1);
normals = zeros(vertices,2);
tangents = zeros(vertices,2);

%while time < maxTime
while count < 1
    count = count + 1;
    
    %errors = zeros(1,size(points,1));
    
    % iterate over each point in the circle
    for ii = 1:size(points,1)
        % find the flow strength
        
        [un,normal,tangent] = calculateIntegral(ii,points);
        uns(ii) = un;
        normals(ii,:) = normal;
        tangents(ii,:) = tangent;
        
%         phiCorrect = atan2(points(ii,2),points(ii,1));
%         normal = [cos(phiCorrect),sin(phiCorrect)];
%         tangent = [-sin(phiCorrect),cos(phiCorrect)];
        
%         news = [points(ii,:);points(ii,:)+0.3*tangent];
%         plot(news(:,1),news(:,2),'r-x')
%         
%         news = [points(ii,:);points(ii,:)+0.3*normal];
%         plot(news(:,1),news(:,2),'k-x')
        
        
        % this has to be applied in the direction of the outward pointing
        % normal
        
        flow = un.*normal;

        % move that point by how much
        
        newPoints(ii,:) = points(ii,:) + flow*dt;
        
        errors(ii) = findDist(points(ii,:),newPoints(ii,:));
        
        
%         if ii == size(points,1)-1 || ii == size(points,1)-2
%             disp('old')
%             points(ii,:)
%             disp('new')
%             newPoints(ii,:)
%             disp('nonsense')
%             un
%             normal
% 
%         end
           
    end
    
    % the error is just un actually
    %error = mean(errors);
    
    points = newPoints;
    
    figure;plot(angles,errors)
    title('errors')
    
    figure;plot(angles,uns)
    title('uns')
    
    max(uns)-min(uns)
    
    figure; hold on;
    plot(angles,abs(cos(angles)-normals(:,1)),'r-x')
%     plot(angles,cos(angles),'r-x');
%     plot(angles,normals(:,1),'k-x')
    title('x component theta hat')
    
    figure; hold on;
    plot(angles,abs(sin(angles)-normals(:,2)),'r-x')
%     plot(angles,sin(angles),'r-x');
%     plot(angles,normals(:,2),'k-x')
    title('y component theta hat')
    
    figure; hold on;
    plot(angles,abs(-sin(angles)-tangents(:,1)),'r-x')
%     plot(angles,-sin(angles),'r-x');
%     plot(angles,tangents(:,1),'k-x')
    title('x component r hat')
    
    figure; hold on;
    plot(angles,abs(cos(angles)-tangents(:,2)),'r-x')
%     plot(angles,cos(angles),'r-x');
%     plot(angles,tangents(:,2),'k-x')
    title('y component r hat')

    % show us the picture        
%     if mod(count,3) == 0
%         fig = figure;
%         plot(points(:,1),points(:,2),'x');
%         xlim([-1.2,1.2]);ylim([-2.2,2.2]);
%         
%         frame = getframe(fig);
%         im{idx} = frame2im(frame);
%         idx = idx + 1;
%     end
    %pause();
    
    time = time + dt;

end

% saveLocation = 'testCircle.tif';
% imwrite(im{1},saveLocation)
% for i=2:numel(im)
%     imwrite(im{i},saveLocation,'WriteMode','append')
% end

end

function [result,normal,properTangent] = calculateIntegral(index,points)
% calculate the integral at the point points(index) = r0

% calculate the normal at r0

r0 = points(index,:);

% add the end of points to the beginning and beginning to the end so it
% loops round (makes it more easy to do the loop later)
% it doesn't matter where we start the integral. And repeats of r0 are
% annoying, so arrange it so this doesn't happen (p124 lab book)

% we want to arrange it so that our r0 is at the beginning/end of the
% circle


if index > 1 && index < size(points,1)
    % we are not at the end points
    points2 = [points(index-1:end,:);points(1:index+1,:)];
elseif index == 1
    % at the start
    points2 = [points(end,:);points;points(1:2,:)];
else
    % at the end
    points2 = [points(end-1:end,:);points;points(1,:)];
end

% our new index is always 2

newIndex = 2;

[normal,properTangent] = findTangentQuadratic(points2(newIndex-1,:),points2(newIndex,:),points2(newIndex+1,:),1);
% 
% phiCorrect = atan2(points2(newIndex,2),points2(newIndex,1));
% normal = [cos(phiCorrect),sin(phiCorrect)];
% properTangent = [-sin(phiCorrect),cos(phiCorrect)];

integrandValsX = zeros(size(points2,1)-2,1);
integrandValsY = zeros(size(points2,1)-2,1);

distances = zeros(size(points2,1)-2,1);

% note we loop from 2 to size - 1 here as we added those three extra points
% on the ends of our original points vector
for i=2:(size(points2,1)-1)
    
    % calculate the tangent at points(i)
    [~,tangent] = findTangentQuadratic(points2(i-1,:),points2(i,:),points2(i+1,:),1);
    
%     phi = atan2(points2(i,2),points2(i,1));
%     tangent = [-sin(phi),cos(phi)];

    
    %  find the value of the integrand at points(i)
    if i == 2 || i==size(points2,1)-1
        integrandX = 0;
        integrandY = 0;
    else
        [integrandX,integrandY] = calcIntegrand2(r0,points2(i,:),tangent);
    end
    
    integrandValsX(i-1) = integrandX/(2*pi);
    integrandValsY(i-1) = integrandY/(2*pi);
    
    % now that we have the values for the integral we need to know the
    % distances between points. These get added together so we know the
    % total distance around
    
    if i==2 
        % if we are at our starting point set it to be zero
        distances(i-1) = 0;
    else
        distances(i-1) = distances(i-2) + findDist(points2(i,:),points2(i-1,:));
    end
    
end

% do the integration
% resultX = trapz(distances,integrandValsX);
% resultY = trapz(distances,integrandValsY);

resultX = irregularSimpson(distances,integrandValsX);
resultY = irregularSimpson(distances,integrandValsY);

if isnan(resultX) || isnan(resultY)
    index
    r0 = points2(index,:)
    
    integrandValsX
    integrandValsY
    
    disp('sum')
    sum(integrandValsX(2:end-1))
    
    
    figure;
    plot(distances,integrandValsX)
    
    error('Was NaN')
    
end


% we have to dot with the normal to get the size of the flow in the normal
% direction
result = resultX*normal(1) + resultY*normal(2);

end

function [integrandX,integrandY] = calcIntegrand(r0,rs,t)
% point calculating tension at, point integrating around, normal at r0,
% tangent at rs

dr = rs - r0;

square = (dr(1).^2 + dr(2).^2);

integrandX = t(1)*t(1)*calcJxxx(dr(1),square);

integrandX = integrandX + t(1)*t(2)*calcSumJxxyJxyx(dr(1),dr(2),square);

integrandX = integrandX + t(2)*t(2)*calcJxyy(dr(1),dr(2),square);

integrandY = t(1)*t(1)*calcJyxx(dr(1),dr(2),square);

integrandY = integrandY + t(1)*t(2)*calcSumJyxyJyyx(dr(1),dr(2),square);

integrandY = integrandY + t(2)*t(2)*calcJyyy(dr(2),square);

end

function [integrandX,integrandY] = calcIntegrand2(r0,rs,t)

dr = rs - r0;

square = (dr(1).^2 + dr(2).^2);

Jexpression = calcJexpression(dr(1),dr(2),t(1),t(2));

integrandX = (dr(1)/(square.^2))*Jexpression;

integrandY = (dr(2)/(square.^2))*Jexpression;


end

function result = calcJexpression(dX,dY,tX,tY)

result = tX*tX*(dY*dY-dX*dX) - 4*tX*tY*dX*dY + tY*tY*(dX*dX-dY*dY);

end

function result = calcJxxx(dX,square)
% dY is not needed to  calculate Jxxx

result = dX/square - 2*(dX.^3)/square.^2;

end

function result = calcSumJxxyJxyx(dX,dY,square)
% We only need the sum of these bits not the bits individually so calculate
% the sum instead

result = - 4*dY*(dX.^2)/square.^2;

end

function result = calcJxyy(dX,dY,square)

result = dX/square - 2*dX*(dY.^2)/square.^2;

end

function result = calcJyxx(dX,dY,square)

result = dY/square - 2*dY*(dX.^2)/square.^2;

end

function result = calcSumJyxyJyyx(dX,dY,square)

result = -4*(dX*dY.^2)/square.^2;

end


function result = calcJyyy(dY,square)
% dX is not needed to  calculate Jyyy

result = dY/square - 2*(dY.^3)/square.^2;

end
