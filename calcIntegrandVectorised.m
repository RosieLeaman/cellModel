function integrands = calcIntegrandVectorised(t,r0,r0rotation,splineX,splineY,plotYes,normalYes,a,b)

% we want to return the integrand at every point p(t)

% first we find the points

points = zeros(numel(t),2);

points(:,1) = ppval(splineX,t);
points(:,2) = ppval(splineY,t);

% we then need to rotate them all to get the correct dx, dy
% unfortunately our points matrix has points as rows, not column vectors,
% so we have to transpose first and then transpose back to keep consistency
points = (r0rotation*points')';

% we must also rotate r0
r0 = (r0rotation*r0')';

% we also need the tangents at these points

[tangents,~] = findTangentFromSplines(t,splineX,splineY,0);

unitTangents = tangents;
for i=1:size(unitTangents,1)
    unitTangents(i,:) = unitTangents(i,:)./norm(unitTangents(i,:));
end

% these also have to be rotated the same as the points
% it may be a better idea to re-find the tangents from the new points, not
% sure atm
tangents = (r0rotation*tangents')';

% then we can work out the integrand easy
integrands = zeros(1,numel(t));
integrandsX = zeros(1,numel(t));
integrandsY = zeros(1,numel(t));
Jexpressions = zeros(1,numel(t));
ds = zeros(1,numel(t));

for i=1:numel(integrands)
    [integrandX,Jexpression] = calcIntegrand(r0,points(i,:),unitTangents(i,:));
    
   % this needs reordering, what on earth is this nonsense atm
    
    % save these, remembering to multiply the integrand by |ds|
    absDs = sqrt(tangents(i,1)^2+tangents(i,2)^2);
    ds(i) = absDs;
    
    %absDs = 1;
    
    integrandsX(i) = integrandX*absDs;
   
    Jexpressions(i) = Jexpression;
    
    integrands(i) = integrandsX(i);
end

% when parametrised by theta, for an ellipse with x stretch a and y stretch
% b, the correct comparison functions are given by:

% note that these are parametrised so that the integral range is 0-1 NOT
% 0-2pi

% ENTIRE INTEGRAL:
% fun5 = @(z) (2*pi*a*(a^2*sin(2*pi*z).^2 + b^2*cos(2*pi*z).^2).^(1/2).*(cos(2*pi*z) - 1).*((a^2.*sin(2*pi*z).^2.*(b^2*sin(2*pi*z).^2 - a^2.*(cos(2*pi*z) - 1).^2))./(a^2.*sin(2*pi*z).^2 + b^2.*cos(2*pi*z).^2) - (b^2*cos(2*pi*z).^2.*(b^2.*sin(2*pi*z).^2 - a^2.*(cos(2*pi*z) - 1).^2))./(a^2.*sin(2*pi*z).^2 + b^2.*cos(2*pi*z).^2) + (4*a^2*b^2.*cos(2*pi*z).*sin(2*pi*z).^2.*(cos(2*pi*z) - 1))./(a^2.*sin(2*pi*z).^2 + b^2.*cos(2*pi*z).^2)))./(b^2.*sin(2*pi*z).^2 + a^2.*(cos(2*pi*z) - 1).^2).^2;

% x and y coordinates
% xvalues = @(z) a*cos(2*pi*z);
% yvalues = @(z) b*sin(2*pi*z);

% tangent and unit tangents
% tangentXunit = @(z) -a*sin(2*pi*z)./sqrt((a*sin(2*pi*z)).^2+(b*cos(2*pi*z)).^2);
% tangentYunit = @(z) b*cos(2*pi*z)./sqrt((a*sin(2*pi*z)).^2+(b*cos(2*pi*z)).^2);
% 
% tangentX = @(z) -2*pi*a*sin(2*pi*z);
% tangentY = @(z) 2*pi*b*cos(2*pi*z);

% ds = sqrt(dx/dt^2+dy/dt^2)
% funDs = @(z) 2*pi*(a^2*sin(2*pi*z).^2 + b^2*cos(2*pi*z).^2).^(1/2);


if plotYes==1
    figure; hold on;
    plot(points(:,1),points(:,2),'x-')
end

if normalYes == 1
    disp('x integrand only')
    integrands = integrandsX;
elseif normalYes == 2
    disp('y integrand only')
    integrands = integrandsY;
end

end

function [integrandX,Jexpression] = calcIntegrand(r0,rs,t)
% in the case where we rotate, we only need to calculate integrandX 

dr = rs - r0;

square = (dr(1).^2 + dr(2).^2);

Jexpression = calcJexpression(dr(1),dr(2),t(1),t(2));

integrandX = (dr(1)/(square.^2))*Jexpression;

%integrandY = (dr(2)/(square.^2))*Jexpression;

end

function result = calcJexpression(dX,dY,tX,tY)

result = tX*tX*(dY*dY-dX*dX) - 4*tX*tY*dX*dY + tY*tY*(dX*dX-dY*dY);

end

function [integrandX,integrandY,Jexpression,checkX,checkY] = calcIntegrand2(r0,rs,t)

dr = rs - r0;

square = (dr(1).^2 + dr(2).^2);

Jexpression = calcJexpression2(dr(1),dr(2),t(1),t(2),square);

integrandX = (dr(1)/(square))*Jexpression;

integrandY = (dr(2)/(square))*Jexpression;

checkX = dr(1)/(square.^2);
checkY = dr(2)/(square.^2);

end

function result = calcJexpression2(dX,dY,tX,tY,square)

result = tX*tX*(1-2*(dX*dX/square)) - 4*tX*tY*(dX*dY/square) + tY*tY*(1-2*(dY*dY/square));

end
