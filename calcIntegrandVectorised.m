function integrands = calcIntegrandVectorised(t,r0,r0normal,splineX,splineY,plotYes,normalYes)

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

% analytic tangent for a circle
% tangents(:,1) = -points(:,2);
% tangents(:,2) = points(:,1);

% make it a unit tangent
unitTangents = tangents;
for i=1:size(unitTangents,1)
    unitTangents(i,:) = unitTangents(i,:)/norm(unitTangents(i,:));
end

% then we can work out the integrand easy
integrands = zeros(1,numel(t));
integrandsX = zeros(1,numel(t));
integrandsY = zeros(1,numel(t));
Jexpressions = zeros(1,numel(t));
checksX = zeros(1,numel(t));
checksY = zeros(1,numel(t));
ds = zeros(1,numel(t));

for i=1:numel(integrands)
    [integrandX,integrandY,Jexpression,checkX,checkY] = calcIntegrand(r0,points(i,:),unitTangents(i,:));
    
    integrands(i) = integrandX*r0normal(1) + integrandY*r0normal(2);
    
    % save these, remembering to multiply the integrand by |ds|
    absDs = sqrt(tangents(i,1)^2+tangents(i,2)^2);
    ds(i) = absDs;
    %absDs = 1;
    integrandsX(i) = integrandX*absDs;
    integrandsY(i) = integrandY*absDs;
    Jexpressions(i) = Jexpression;
    checksX(i) = checkX;
    checksY(i) = checkY;
end

if plotYes==1
    figure;
    plot(points(:,1),points(:,2),'x-')
    title(['points, r0 is (',num2str(r0(1)),',',num2str(r0(2)),')'])
%     
%     figure; hold on;
%     plot(points(:,1),points(:,2),'x-')
%     for i=1:1000:size(tangents,1)
%         plot([points(i,1),points(i,1)+0.5*tangents(i,1)],[points(i,2),points(i,2)+0.5*tangents(i,2)],'ro-')
%     end
    
%     figure;
%     subplot(1,3,1)
%     plot(t,integrands,'x-')
%     title('integrand')
%     subplot(1,3,2)
%     plot(points(:,1),integrandsX,'x-')
%     title('integrand X')
%     subplot(1,3,3)
%     plot(points(:,2),integrandsY,'x-')
%     title('integrand Y')
%     
%     figure;
%     subplot(1,3,1)
%     plot(points(:,1),Jexpressions,'x')
%     title('J')
%     subplot(1,3,2)
%     plot(points(:,1),checksX,'x')
%     title('x prefactor')
%     subplot(1,3,3)
%     plot(points(:,2),checksY,'x')
%     title('y prefactor')
    
    a = 2;
    b = 2;
    fun5 = @(z) (a*(a^2*sin(z).^2 + b^2*cos(z).^2).^(1/2).*(cos(z) - 1).*((a^2.*sin(z).^2.*(b^2*sin(z).^2 - a^2.*(cos(z) - 1).^2))./(a^2.*sin(z).^2 + b^2.*cos(z).^2) - (b^2*cos(z).^2.*(b^2.*sin(z).^2 - a^2.*(cos(z) - 1).^2))./(a^2.*sin(z).^2 + b^2.*cos(z).^2) + (4*a^2*b^2.*cos(z).*sin(z).^2.*(cos(z) - 1))./(a^2.*sin(z).^2 + b^2.*cos(z).^2)))./(b^2.*sin(z).^2 + a^2.*(cos(z) - 1).^2).^2;
    
    figure;hold on;
    plot(t,integrandsX)
    plot(t,fun5(t))
    title('whole integrand')
    
    % check the tangents
    
    xvalues = @(z) a*cos(z);
    
    tangentX = @(z) -a*sin(z);
    tangentY = @(z) b*cos(z);

    figure;hold on;
    plot(t,tangents(:,1))
    plot(t,tangentX(t))
    title('tangent x')
    
    figure;hold on;
    plot(t,tangents(:,2))
    plot(t,tangentY(t))
    title('tangent y')
    
    figure;
    subplot(2,2,1)
    plot(t,xvalues(t))
    title('correct x values')
    subplot(2,2,2)
    plot(t,tangentX(t))
    title('correct x derivatives')
    subplot(2,2,3)
    plot(t,points(:,1))
    title('x values')
    subplot(2,2,4)
    plot(t,ppval(splineXdiff,t))
    title('spline derivative')
    
    funDs = @(z) (a^2*sin(z).^2 + b^2*cos(z).^2).^(1/2);
    
    figure; hold on;
    plot(t,ds)
    plot(t,funDs(t))
    title('diffs ds')
    
    funJ = @(z) (a^2.*sin(z).^2.*(b^2*sin(z).^2 - a^2.*(cos(z) - 1).^2))./(a^2.*sin(z).^2 + b^2.*cos(z).^2) - (b^2*cos(z).^2.*(b^2.*sin(z).^2 - a^2.*(cos(z) - 1).^2))./(a^2.*sin(z).^2 + b^2.*cos(z).^2) + (4*a^2*b^2.*cos(z).*sin(z).^2.*(cos(z) - 1))./(a^2.*sin(z).^2 + b^2.*cos(z).^2);

    figure;hold on;
    plot(t,Jexpressions)
    plot(t,funJ(t))
    title('diffs J')
    
    figure; hold on;
    plot(points(:,1)-r0(1))
    plot(points(:,2)-r0(2))
    title('dr')
    
    funPrefactorX = @(z) (a*(cos(z) - 1))./(b^2.*sin(z).^2 + a^2.*(cos(z) - 1).^2).^2;
    
    figure;hold on;
    plot(t,checksX)
    plot(t,funPrefactorX(t))
    title('diffs prefactor')
    
    funJplusPre = @(z) (a*(cos(z) - 1).*((a^2.*sin(z).^2.*(b^2*sin(z).^2 - a^2.*(cos(z) - 1).^2))./(a^2.*sin(z).^2 + b^2.*cos(z).^2) - (b^2*cos(z).^2.*(b^2.*sin(z).^2 - a^2.*(cos(z) - 1).^2))./(a^2.*sin(z).^2 + b^2.*cos(z).^2) + (4*a^2*b^2.*cos(z).*sin(z).^2.*(cos(z) - 1))./(a^2.*sin(z).^2 + b^2.*cos(z).^2)))./(b^2.*sin(z).^2 + a^2.*(cos(z) - 1).^2).^2;
    
    figure; hold on;
    plot(t,checksX.*Jexpressions)
    plot(t,funJplusPre(t))
    title('diffs jxpre')
    
end

if normalYes == 1
    integrands = integrandsX;
elseif normalYes == 2
    integrands = integrandsY;
end

end

function [integrandX,integrandY,Jexpression,checkX,checkY] = calcIntegrand(r0,rs,t)

dr = rs - r0;

square = (dr(1).^2 + dr(2).^2);

Jexpression = calcJexpression(dr(1),dr(2),t(1),t(2));

integrandX = (dr(1)/(square.^2))*Jexpression;

integrandY = (dr(2)/(square.^2))*Jexpression;

checkX = dr(1)/(square.^2);
checkY = dr(2)/(square.^2);

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

