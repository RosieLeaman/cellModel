function surfaceTensionTest()

% get a circle
points = findVerticesNewMaterial([0,0],101,1);

% for a certain length of time

time = 0;
maxTime = 0.05;
dt = 0.01;

figure;
plot(points(:,1),points(:,2),'x-');
axis square
hold on;

pause();

while time < maxTime
    
    newPoints = points;

    % iterate over each point in the circle
    for ii = 1:size(points,1)
        % find the flow strength
        
        [un,normal] = calculateIntegral(ii,points);
        
        % this has to be applied in the direction of the outward pointing
        % normal
        
        if ii== 26 || ii== 76
            time
            ii
            un
            normal            
        end
        
        flow = un.*normal;

        % move that point by how much
        
        newPoints(ii,:) = points(ii,:) + flow*dt;
           
    end
    
    points = newPoints;

    % show us the picture
    
    figure;
    plot(points(:,1),points(:,2),'x');
    axis square
    
    %pause();
    
    time = time + dt;

end

end

function [result,normal] = calculateIntegral(index,points)
% calculate the integral at the point points(index) = r0

% calculate the normal at r0 which is used every time

r0 = points(index,:);

% add the end of points to the beginning and beginning to the end so it
% loops round (makes it more easy to do the loop later)

points2 = [points(end,:);points;points(1,:)];

% this shifts our original index up by 1
newIndex = index + 1;

[normal,~] = findTangentQuadratic(points2(newIndex-1,:),points2(newIndex,:),points2(newIndex+1,:),1);

integrandVals = zeros(size(points2,1)-2,1);

distances = zeros(size(points2,1)-2,1);


% note we loop from 2 to size - 1 here as we added those two extra points
% on the ends of our original points vector
for i=2:(size(points2,1)-1)
    
    % calculate the tangent at points(i)
    [~,tangent] = findTangentQuadratic(points2(i-1,:),points2(i,:),points2(i+1,:),1);

    
    %  find the value of the integrand at points(i)
    if i ~= newIndex
        integrand = calcIntegrand(r0,points2(i,:),normal,tangent);
    else
        integrand = 0;
    end
    
    integrandVals(i-1) = integrand;
    
    % now that we have the values for the integral we need to know the
    % distances between points. These get added together so we know the
    % total distance around
    
    if i==2
        distances(i-1) = 0;
    else
        distances(i-1) = distances(i-2) + findDist(points2(i,:),points2(i-1,:));
    end
end

result = trapz(distances,integrandVals);

end

function integrand = calcIntegrand(r0,rs,n,t)
% point calculating tension at, point integrating around, normal at r0,
% tangent at rs

dr = rs - r0;

square = (dr(1).^2 + dr(2).^2);

integrand = n(1)*t(1)*t(1)*calcJxxx(dr(1),dr(2));

integrand = integrand + (2*n(1)*t(1)*t(2) + n(2)*t(1)*t(1))*calcJyxx(dr(1),dr(2));

integrand = integrand + (2*n(2)*t(1)*t(2) + n(1)*t(2)*t(2))*calcJxyy(dr(1),dr(2));

integrand = integrand + n(2)*t(2)*t(2)*calcJyyy(dr(1),dr(2));

end

function integrand = calcIntegrandY(r0,rs,n,t)
% point calculating tension at, point integrating around, normal at r0,
% tangent at rs

dr = rs - r0;

square = (dr(1).^2 + dr(2).^2);

integrand = n(1)*t(1)*t(1)*calcJxxx(dr(1),dr(2));

integrand = integrand + (2*n(1)*t(1)*t(2) + n(2)*t(1)*t(1))*calcJyxx(dr(1),dr(2));

integrand = integrand + (2*n(2)*t(1)*t(2) + n(1)*t(2)*t(2))*calcJxyy(dr(1),dr(2));

integrand = integrand + n(2)*t(2)*t(2)*calcJyyy(dr(1),dr(2));

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

result = dY/square - 2*(dY*dX.^2)/square.^2;

end

function result = calcSumJyxyJyyx(dX,dY,square)

result = -4*(dX*dY.^2)/square.^2;

end


function result = calcJyyy(dY,square)
% dX is not needed to  calculate Jyyy

result = dY/square - 2*(dY.^3)/square.^2;

end