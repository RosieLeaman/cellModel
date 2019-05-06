function integrands = calcIntegrandVectorised(t,r0,r0normal,splineX,splineY,plotYes)

% we want to return the integrand at every point p(t)

% first we find the points

points = zeros(numel(t),2);

points(:,1) = ppval(splineX,t);
points(:,2) = ppval(splineY,t);

% we also need the tangents at these points

splineXdiff = fnder(splineX);
splineYdiff = fnder(splineY);

xdiffs = ppval(splineXdiff,t);
ydiffs = ppval(splineYdiff,t);

tangents = zeros(numel(xdiffs),2);
tangents(:,1) = xdiffs;
tangents(:,2) = ydiffs;

% make it a unit tangent
for i=1:size(tangents,1)
    tangents(i,:) = tangents(i,:)/norm(tangents(i,:));
end


% then we can work out the integrand easy

integrands = zeros(1,numel(t));
integrandsX = zeros(1,numel(t));
integrandsY = zeros(1,numel(t));

for i=1:numel(integrands)
    [integrandX,integrandY] = calcIntegrand(r0,points(i,:),tangents(i,:),plotYes);
    
    integrands(i) = integrandX*r0normal(1) + integrandY*r0normal(2);
    integrandsX(i) = integrandX;
    integrandsY(i) = integrandY;
end

if plotYes==1
    integrandsX
    integrandsY
    points(1,:)
    points(3,:)-r0
    figure;
    plot(points(:,1),points(:,2),'x')
    figure;
    subplot(1,2,1)
    plot(points(:,1),integrandsX,'x-')
    subplot(1,2,2)
    plot(points(:,2),integrandsY,'x-')
end

end

function [integrandX,integrandY] = calcIntegrand(r0,rs,t,plotYes)

dr = rs - r0;

square = (dr(1).^2 + dr(2).^2);

Jexpression = calcJexpression(dr(1),dr(2),t(1),t(2));

integrandX = (dr(1)/(square.^2))*Jexpression;

integrandY = (dr(2)/(square.^2))*Jexpression;

end

function result = calcJexpression(dX,dY,tX,tY)

result = tX*tX*(dY*dY-dX*dX) - 4*tX*tY*dX*dY + tY*tY*(dX*dX-dY*dY);

end
