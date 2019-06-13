function integrands = calcIntegrandVectorised(t,r0,r0normal,splineX,splineY,plotYes,normalYes,a,b)

% we want to return the integrand at every point p(t)

% first we find the points

points = zeros(numel(t),2);

points(:,1) = ppval(splineX,t);
points(:,2) = ppval(splineY,t);

% we also need the tangents at these points

[tangents,~] = findTangentFromSplines(t,splineX,splineY,0);

unitTangents = tangents;
for i=1:size(unitTangents,1)
    unitTangents(i,:) = unitTangents(i,:)./norm(unitTangents(i,:));
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
    
    fun5 = @(z) (2*pi*a*(a^2*sin(2*pi*z).^2 + b^2*cos(2*pi*z).^2).^(1/2).*(cos(2*pi*z) - 1).*((a^2.*sin(2*pi*z).^2.*(b^2*sin(2*pi*z).^2 - a^2.*(cos(2*pi*z) - 1).^2))./(a^2.*sin(2*pi*z).^2 + b^2.*cos(2*pi*z).^2) - (b^2*cos(2*pi*z).^2.*(b^2.*sin(2*pi*z).^2 - a^2.*(cos(2*pi*z) - 1).^2))./(a^2.*sin(2*pi*z).^2 + b^2.*cos(2*pi*z).^2) + (4*a^2*b^2.*cos(2*pi*z).*sin(2*pi*z).^2.*(cos(2*pi*z) - 1))./(a^2.*sin(2*pi*z).^2 + b^2.*cos(2*pi*z).^2)))./(b^2.*sin(2*pi*z).^2 + a^2.*(cos(2*pi*z) - 1).^2).^2;

    figure;hold on;
    plot(t,integrandsX,'x')
    plot(t,fun5(t),'o')
    title('whole integrand')
    
    diffs = zeros(1,numel(integrandsX));
    for i=1:numel(diffs)
        diffs(i) = integrandsX(i)-fun5(t(i));
    end
    figure;
    plot(t,diffs);
    title('integrand diffs')
    
    % check the tangents
    
    xvalues = @(z) a*cos(2*pi*z);
    yvalues = @(z) b*sin(2*pi*z);
    
    tangentXunit = @(z) -a*sin(2*pi*z)./sqrt((a*sin(2*pi*z)).^2+(b*cos(2*pi*z)).^2);
    tangentYunit = @(z) b*cos(2*pi*z)./sqrt((a*sin(2*pi*z)).^2+(b*cos(2*pi*z)).^2);
    
    tangentX = @(z) -2*pi*a*sin(2*pi*z);
    tangentY = @(z) 2*pi*b*cos(2*pi*z);
    
    figure;hold on;
    plot(t,tangents(:,1),'x')
    plot(t,tangentX(t),'o')
    title('tangent x')
    
    figure;hold on;
    plot(t,tangents(:,2),'x')
    plot(t,tangentY(t),'o')
    title('tangent y')

    figure;hold on;
    plot(t,unitTangents(:,1),'x')
    plot(t,tangentXunit(t),'o')
    title('unit tangent x')
    
    figure;hold on;
    plot(t,unitTangents(:,2),'x')
    plot(t,tangentYunit(t),'o')
    title('unit tangent y')
    
    funDs = @(z) (a^2*sin(2*pi*z).^2 + b^2*cos(2*pi*z).^2).^(1/2);
    
    diffs = zeros(1,numel(ds));
    for i=1:numel(diffs)
        diffs(i) = ds(i)-funDs(t(i));
    end
    figure;
    plot(t,diffs);
    title('ds diffs')
%     
%     funJ = @(z) (a^2.*sin(z).^2.*(b^2*sin(z).^2 - a^2.*(cos(z) - 1).^2))./(a^2.*sin(z).^2 + b^2.*cos(z).^2) - (b^2*cos(z).^2.*(b^2.*sin(z).^2 - a^2.*(cos(z) - 1).^2))./(a^2.*sin(z).^2 + b^2.*cos(z).^2) + (4*a^2*b^2.*cos(z).*sin(z).^2.*(cos(z) - 1))./(a^2.*sin(z).^2 + b^2.*cos(z).^2);
% 
%     figure;hold on;
%     plot(t,Jexpressions,'x')
%     plot(t,funJ(2*pi*t),'o')
%     title('diffs J')
%     
%     diffs = zeros(1,numel(Jexpressions));
%     for i=1:numel(diffs)
%         diffs(i) = Jexpressions(i)-funJ(2*pi*t(i));
%     end
%     figure;
%     plot(t,diffs);
%     title('J diffs')
%     
% %     figure; hold on;
% %     plot(t,points(:,1)-r0(1))
% %     plot(t,points(:,2)-r0(2))
% %     title('dr')
%     
%     funPrefactorX = @(z) (a*(cos(z) - 1))./(b^2.*sin(z).^2 + a^2.*(cos(z) - 1).^2).^2;
%     
%     figure;hold on;
%     plot(t,checksX,'x')
%     plot(t,funPrefactorX(2*pi*t),'o')
%     title('diffs prefactor')
%     
%     diffs = zeros(1,numel(checksX));
%     for i=1:numel(diffs)
%         diffs(i) = checksX(i)-funPrefactorX(2*pi*t(i));
%     end
%     figure;
%     plot(t,diffs);
%     title('prefactor diffs')
%     
%     funJplusPre = @(z) (a*(cos(z) - 1).*((a^2.*sin(z).^2.*(b^2*sin(z).^2 - a^2.*(cos(z) - 1).^2))./(a^2.*sin(z).^2 + b^2.*cos(z).^2) - (b^2*cos(z).^2.*(b^2.*sin(z).^2 - a^2.*(cos(z) - 1).^2))./(a^2.*sin(z).^2 + b^2.*cos(z).^2) + (4*a^2*b^2.*cos(z).*sin(z).^2.*(cos(z) - 1))./(a^2.*sin(z).^2 + b^2.*cos(z).^2)))./(b^2.*sin(z).^2 + a^2.*(cos(z) - 1).^2).^2;
%     
%     figure; hold on;
%     plot(t,checksX.*Jexpressions)
%     plot(t,funJplusPre(2*pi*t))
%     title('diffs jxpre')
    
end

if normalYes == 1
    disp('x integrand only')
    integrands = integrandsX;
elseif normalYes == 2
    disp('y integrand only')
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

