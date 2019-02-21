function surfaceTensionTest()

% get a circle
vertices = 251;
radius = 1;
amount = (vertices*radius/2)*sin(2*pi/vertices);
points = findVerticesNewMaterial([0,0],101,amount);
newPoints = points;

% for a certain length of time

time = 0;
maxTime = 0.03;
dt = 0.01;

figure;
plot(points(:,1),points(:,2),'x-');
axis square

%pause();

count = 0;
while time < maxTime
    count = count + 1;
    
    % iterate over each point in the circle
    for ii = 1:size(points,1)
        % find the flow strength
        
        [un,normal] = calculateIntegral(ii,points);
        
        % this has to be applied in the direction of the outward pointing
        % normal
        
        flow = un.*normal;

        % move that point by how much
        
        newPoints(ii,:) = points(ii,:) + flow*dt;
        
        if ii == 1
            disp('old')
            points(ii,:)
            disp('new')
            newPoints(ii,:)
            disp('nonsense')
            un
            normal
            
            phi = atan2(points(ii,2),points(ii,1))
            trueNormal = [cos(phi),sin(phi)]
            trueTangent = [-sin(phi),cos(phi)]
        end
           
    end
    
    points = newPoints;

    % show us the picture        
    if mod(count,3) == 0
        figure;
        plot(points(:,1),points(:,2),'x');
        axis square
    end
    %pause();
    
    time = time + dt;
    
    points

end

end

function [result,normal] = calculateIntegral(index,points)
% calculate the integral at the point points(index) = r0

% calculate the normal at r0

r0 = points(index,:);

% add the end of points to the beginning and beginning to the end so it
% loops round (makes it more easy to do the loop later)
% it doesn't matter where we start the integral. And we can't allow repeats

points2 = [points(end,:);points;points(1,:);points(2,:)];

% this shifts our original index up by 1
newIndex = index + 1;

[normal,~] = findTangentQuadratic(points2(newIndex-1,:),points2(newIndex,:),points2(newIndex+1,:),1);

%phi = atan2(points2(newIndex,2),points2(newIndex,1));
%normal = [cos(phi),sin(phi)];

integrandValsX = zeros(size(points2,1)-2,1);
integrandValsY = zeros(size(points2,1)-2,1);

distances = zeros(size(points2,1)-2,1);

% note we loop from 2 to size - 1 here as we added those two extra points
% on the ends of our original points vector
for i=2:(size(points2,1)-1)
    
    % calculate the tangent at points(i)
    [~,tangent] = findTangentQuadratic(points2(i-1,:),points2(i,:),points2(i+1,:),1);
    
    %phi = atan2(points2(i,2),points2(i,1));
    %tangent = [-sin(phi),cos(phi)];

    
    %  find the value of the integrand at points(i)
    if i ~= newIndex
        [integrandX,integrandY] = calcIntegrand2(r0,points2(i,:),tangent);
        if index == 1 && i == (size(points2,1)-1)
            disp('ehh')
            integrandX
            integrandY
            r0
            points2(i,:)
            tangent
            disp('ehhhhhh')
        end
    else
        integrandX = 0;
        integrandY = 0;
    end
    
    integrandValsX(i-1) = integrandX;
    integrandValsY(i-1) = integrandY;
    
    % now that we have the values for the integral we need to know the
    % distances between points. These get added together so we know the
    % total distance around
    
    if i==2
        distances(i-1) = 0;
    else
        distances(i-1) = distances(i-2) + findDist(points2(i,:),points2(i-1,:));
    end
end

% do the integration
resultX = trapz(distances,integrandValsX);
resultY = trapz(distances,integrandValsY);

if index == 1
    integrandValsX
    integrandValsY
    distances
    resultX
    resultY
    dr = points2(1,:) - r0

    square = (dr(1).^2 + dr(2).^2)    
    
    index2 = (size(points2,1)-1);
    [~,tangent] = findTangentQuadratic(points2(index2-1,:),points2(index2,:),points2(index2+1,:),1)
    
    phi = atan2(points2(index2,2),points2(index2,1))
    trueNormal = [cos(phi),sin(phi)]
    trueTangent = [-sin(phi),cos(phi)]
    
    calcJexpression(dr(1),dr(2),tangent(1),tangent(2))
    calcJexpression(dr(1),dr(2),trueTangent(1),trueTangent(2))
    
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