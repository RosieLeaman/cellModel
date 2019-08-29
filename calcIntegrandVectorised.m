% this function calculates the surface tension integrand at the given point
% r0. It expects a vector of parameter values, t
% the parametrisation is EXPECTED to be BY ARC LENGTH

% INPUTS:
% t; 1xn vector of parameter values at which to calculate the integrand. I
%    believe that it is expected to be a row vector
% r0; 1x2 vector. the specific point [r0x,r0y] (point, not parameter value!)
%    at which the integral is to be calculated. 
% r0rotation; 2x2 rotation matrix which rotates the plane so that the point
%    r0 has unit normal (1,0) and unit tangent (0,1)
% splineX; a spline object which interpolates the x co-ordinates given
%    parameter values
% splineY; a spline object which interpolates the y co-ordinates given
%    parameter values
% plotYes; 0 or 1 - 0 no plotting, 1 plotting

%OUTPUTS

% integrands; a 1xt row vector containing the integrand at each input point
%             in t

function integrands = calcIntegrandVectorised(t,r0,r0rotation,splineX,splineY,plotYes)

% we want to return the integrand at every point p(t)

% first we find the points which we do by using the splines to interpolate
points = zeros(numel(t),2);

points(:,1) = ppval(splineX,t);
points(:,2) = ppval(splineY,t);

if plotYes == 1
    disp('before rotate')
    points(1,:)
end

% we then need to rotate them all to get the correct dx, dy
% unfortunately our points matrix has points as rows, not column vectors,
% so we have to transpose first and then transpose back to keep consistency
points = (r0rotation*points')';

if plotYes == 1
    disp('after rotate')
    points(1,:)
end

% we must also rotate r0
r0 = (r0rotation*r0')';

% we also need the tangents at these points which we get from the original
% splines.
[tangents,~] = findTangentFromSplines(t,splineX,splineY,0);

if plotYes == 1
    disp('tangents')
    tangents(1,:)
end

unitTangents = tangents;
for i=1:size(unitTangents,1)
    unitTangents(i,:) = unitTangents(i,:)./norm(unitTangents(i,:));
end

if plotYes == 1
    disp('unit tangents')
    unitTangents(1,:)
end

% these also have to be rotated the same as the points
% it may be a better idea to re-find the tangents from the new points, not
% sure atm
% the non-unit tangents are needed if having to calculate ds
%rotatedTangents = (r0rotation*tangents')';
rotatedUnitTangents = (r0rotation*unitTangents')';

if plotYes == 1
    disp('rotated unit tangents')
    rotatedUnitTangents(1,:)
end

if plotYes == 1
    % plot the tangents to check they look legit
    figure; hold on;
    plot(points(:,1),points(:,2),'x-')
    for i = 1:25:size(points,1)
        plot([points(i,1),points(i,1)+rotatedUnitTangents(i,1)],[points(i,2),points(i,2)+rotatedUnitTangents(i,2)],'o-')
    end
    title('points and tangents')
    figure; hold on;
    plot(points(:,1),points(:,2),'x-')
    title('points')
end

% then we can work out the integrand easy
integrands = zeros(1,numel(t));

for i=1:numel(integrands)
    % note it is OK to have r0 and points here as we rotated them earlier
    integrandX = calcIntegrand(r0,points(i,:),rotatedUnitTangents(i,:));

    % sometimes might need to multiply by ds
    %absDs = sqrt(rotatedTangents(i,1)^2+rotatedTangents(i,2)^2);
    
    % however if we parametrise by arc length s, then ds = 1, so we don't
    % need to
    
    integrands(i) = integrandX/(-2*pi);
    % we have to divide by -2pi here as the J expression calculated below
    % is actually -2pi*J
   
end

if plotYes == 2
    figure;
    plot(integrands)
end

end

function integrandX = calcIntegrand(r0,rs,t)
% in the case where we rotate, we only need to calculate integrandX

dr = rs - r0;

square = (dr(1).^2 + dr(2).^2);

Jexpression = calcJexpression(dr(1),dr(2),t(1),t(2));

integrandX = (dr(1)/(square.^2))*Jexpression;

%integrandY = (dr(2)/(square.^2))*Jexpression;

% if abs(square) <= 10e-16
%     % we have somehow managed to pick the exact same point
%     % return zero
%     integrandX = 0;
% end

end

function result = calcJexpression(dX,dY,tX,tY)

%result = tX*tX*(dY*dY-dX*dX) - 4*tX*tY*dX*dY + tY*tY*(dX*dX-dY*dY);

result = (tX*tX-tY*tY)*(dY*dY-dX*dX) - 4*tX*tY*dX*dY;

%result = (tX - tY)*(tX + tY)*(dY - dX)*(dY + dX) - 4*tX*tY*dX*dY;

%result = (tY*(dX+dY)-tX*(dY-dX))*(tY*(dX-dY)-tX*(dX+dY));

end