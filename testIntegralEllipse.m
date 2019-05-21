function testIntegralEllipse()

a = 2;
b = 2;
p = b/a;

% UNIT TANGENTS, CORRECT [-ASIN(T),BCOS(T)]
% direct substitution
fun5 = @(t) (a*(a^2*sin(t).^2 + b^2*cos(t).^2).^(1/2).*(cos(t) - 1).*((a^2.*sin(t).^2.*(b^2*sin(t).^2 - a^2.*(cos(t) - 1).^2))./(a^2.*sin(t).^2 + b^2.*cos(t).^2) - (b^2*cos(t).^2.*(b^2.*sin(t).^2 - a^2.*(cos(t) - 1).^2))./(a^2.*sin(t).^2 + b^2.*cos(t).^2) + (4*a^2*b^2.*cos(t).*sin(t).^2.*(cos(t) - 1))./(a^2.*sin(t).^2 + b^2.*cos(t).^2)))./(b^2.*sin(t).^2 + a^2.*(cos(t) - 1).^2).^2;

disp('fun5')
integral(fun5,0,2*pi)
disp('values')
fun5(1e-5)
fun5(1e-10)

% simplified
fun6 = @(t) (a.*(a^4*cos(t).^3 - a^4*cos(t).^2 + b^4*cos(t).^2 + b^4*cos(t).^3 + a^4 - a^2*b^2 - a^4.*cos(t) - 2*a^2*b^2*cos(t).^3 + 3*a^2*b^2*cos(t)))./((- a^2*cos(t).^2 + a^2 + b^2*cos(t).^2).^(1/2).*(a^2 + b^2 - a^2*cos(t) + b^2*cos(t)).^2);

disp('fun6')
integral(fun6,0,2*pi)

% simplified combined powers of cos(t)

fun6Alt = @(t) (a.*((a^4+b^4-2*a^2*b^2)*cos(t).^3 + (b^4-a^4)*cos(t).^2 + (3*a^2*b^2-a^4)*cos(t) + a^4-a^2*b^2))./((a^2 + (b^2-a^2)*cos(t).^2).^(1/2).*(a^2 + b^2 + (b^2-a^2)*cos(t)).^2);

disp('fun6Alt')
integral(fun6Alt,0,2*pi)
disp('values')
fun6Alt(1e-5)
fun6Alt(1e-10)

% with substitution b = pa
fun7 = @(t) (3*p^2*cos(t) - cos(t) - cos(t).^2 + cos(t).^3 - 2*p^2*cos(t).^3 + p^4*cos(t).^2 + p^4*cos(t).^3 - p^2 + 1)./((p^2*cos(t).^2 - cos(t).^2 + 1).^(1/2).*(p^2*cos(t) - cos(t) + p^2 + 1).^2);

disp('fun7')
integral(fun7,0,2*pi)

% with substitution b = pa, combined powers of cos(t)
fun8 = @(t) (1-p^2+(3*p^2-1)*cos(t)+(p^4-1)*cos(t).^2+(p^4-2*p^2+1)*cos(t).^3)./((p^2*cos(t).^2 - cos(t).^2 + 1).^(1/2).*(p^2*cos(t) - cos(t) + p^2 + 1).^2);

disp('fun8')
integral(fun8,0,2*pi)

% slight different direct substitution form 
fun9 = @(t) (a*(cos(t) - 1).*((b^2*cos(t).^2 - a^2*sin(t).^2).*(a*(cos(t) - 1) - (b^2*sin(t).^2)./cos(t)) + 4*a*b^2*cos(t).*sin(t).^2))./((a^2*sin(t).^2 + b^2*cos(t).^2).^(1/2).*(a*(cos(t) - 1) + (b^2*sin(t).^2)./cos(t)).^2);

disp('fun9')
integral(fun9,0,2*pi)
disp('values')
fun9(1e-5)
fun9(1e-10)

x = linspace(0,2*pi,1000);

disp('error')
disp(abs(integral(fun6Alt,0,2*pi)-integral(fun5,0,2*pi)));


fun10 = @(x,t) (2*a*sin(t/2 - x/2).*sin(t/2 + x/2).*(a^2*sin(t).^2 + b^2*cos(t).^2).^(1/2).*((b^2*cos(t).^2.*(4*b^2*cos(t/2 + x/2).^2.*sin(t/2 - x/2).^2 - 4*a^2*sin(t/2 - x/2).^2.*sin(t/2 + x/2).^2))./(a^2*sin(t).^2 + b^2*cos(t).^2) - (a^2*sin(t).^2.*(4*b^2*cos(t/2 + x/2).^2.*sin(t/2 - x/2).^2 - 4*a^2*sin(t/2 - x/2).^2.*sin(t/2 + x/2).^2))./(a^2*sin(t).^2 + b^2*cos(t).^2) + (16*a^2*b^2*cos(t/2 + x/2).*sin(t/2 - x/2).^2.*sin(t/2 + x/2).*cos(t).*sin(t))./(a^2*sin(t).^2 + b^2*cos(t).^2)))./(4*a^2*sin(t/2 - x/2).^2.*sin(t/2 + x/2).^2 + 4*b^2*cos(t/2 + x/2).^2.*sin(t/2 - x/2).^2).^2;

figure;hold on;
plot(x,fun5(x))
newFun = @(t) fun10(0,t);
plot(x,newFun(x))

disp('compare fun 10 fun5')
fun5int = integral(fun5,0,2*pi)
newFun = @(t) fun10(0,t);
fun10int = integral(newFun,0,2*pi)
newFun = @(t) fun10(2*pi,t);
fun10int = integral(newFun,0,2*pi)
disp('error')
abs(fun5int-fun10int)

integralAngles = linspace(0,2*pi,50);
integrals = zeros(numel(integralAngles),1);

for i=1:numel(integrals)
    newFun = @(t) fun10(integralAngles(i),t);
    integrals(i) = integral(newFun,0,2*pi);
end

figure;
plot(integralAngles,integrals)
title('integrals')

% figure;
% plot(x,fun5(x))
% xlim([0,2*pi])
% title('fun5')
% figure;
% plot(x,fun8(x))
% xlim([0,2*pi])
% title('fun8')

x2 = linspace(0,0.00001,1000);
figure; hold on;
plot(x2,fun5(x2),'linewidth',2)
plot(x2,fun6Alt(x2),'linewidth',2)
axis tight
grid minor
legend('Direct substitution','Simplified')
xlabel('theta');ylabel('Integrand')

% try with POINTS now
vertices = 250;

points = zeros(vertices,2);
angles = zeros(vertices,1);

index = 0;
while index < vertices
    theta = 0 + (2*pi*index)/vertices;
    points(index+1,1) = a*cos(theta);
    points(index+1,2) = b*sin(theta);
    angles(index+1) = theta;
    index = index + 1;
end

disp('points based integral')
result = calculateIntegral([a,0],a,b)
disp('diff to fun5')
abs(result - integral(fun5,0,2*pi))
disp('diff to fun6Alt')
abs(result - integral(fun6Alt,0,2*pi))


% figure;
% subplot(1,2,1); hold on;
% x3 = linspace(0,2*pi,100); 
% integrands = calcIntegrandVectorX(x3./(2*pi),[a,0],a,b);
% plot(x3,integrands,'x-','linewidth',2)
% plot(x3,fun5(x3),'x--','linewidth',2)
% set(gca,'FontSize',16)
% axis tight
% grid minor
% legend('anonymous function','as for surface tension')
% xlabel('theta');ylabel('Integrand')
% 
% subplot(1,2,2)
% x3 = linspace(0,2*pi,100000); 
% integrands = calcIntegrandVectorX(x3./(2*pi),[a,0],a,b);
% plot(x3,integrands-fun5(x3),'k-','linewidth',2)
% set(gca,'FontSize',16)
% axis tight
% grid minor
% xlabel('theta');ylabel('Difference')

end

function result = calculateIntegral(r0,a,b)
    fun = @(t)calcIntegrandVectorGivenTangentX(t,r0,a,b);

    result = integral(fun,0,1);
end

function integrands = calcIntegrandVectorGivenTangentX(t,r0,a,b)
    % t goes from 0 to 1   
    integrands = zeros(1,numel(t));
    
    % the theta values then go from 0 to 2pi
    theta = (2*pi).*t;
    
%     points(:,1) = a*cos(theta);
%     points(:,2) = b*sin(theta);
    
%     figure;
%     plot(points(:,1),points(:,2),'x-')
    
    tangents(:,1) = -a*sin(theta);
    tangents(:,2) = b*cos(theta);
    
    % make it a unit tangent
    for i=1:size(tangents,1)
        tangents(i,:) = tangents(i,:)/norm(tangents(i,:));
    end
    
    Jexpressions = zeros(1,numel(t));
    
    p=b/a;
    newFun = @(z) (1-p^2+(3*p^2-1)*cos(z)+(p^4-1)*cos(z).^2+(p^4-2*p^2+1)*cos(z).^3)./((p^2*cos(z).^2 - cos(z).^2 + 1).^(1/2).*(p^2*cos(z) - cos(z) + p^2 + 1).^2);
    newFun2 = @(z) (a*(a^2*sin(z).^2 + b^2*cos(z).^2).^(1/2).*(cos(z) - 1).*((a^2.*sin(z).^2.*(b^2*sin(z).^2 - a^2.*(cos(z) - 1).^2))./(a^2.*sin(z).^2 + b^2.*cos(z).^2) - (b^2*cos(z).^2.*(b^2.*sin(z).^2 - a^2.*(cos(z) - 1).^2))./(a^2.*sin(z).^2 + b^2.*cos(z).^2) + (4*a^2*b^2.*cos(z).*sin(z).^2.*(cos(z) - 1))./(a^2.*sin(z).^2 + b^2.*cos(z).^2)))./(b^2.*sin(z).^2 + a^2.*(cos(z) - 1).^2).^2;


    % actually calculate integrand
    for i=1:numel(integrands)
        %[integrandX,~,Jexpression,~,~] = calcIntegrand(r0,points(i,:),tangents(i,:));
        [integrandX,~,Jexpression,~,~] = calcIntegrand2(theta(i),tangents(i,:),a,b);
%         integrandX = calcIntegrand3(theta(i),tangents(i,:),a,b);
%         Jexpression = 0;
%         integrandX = calcIntegrand4(theta(i),tangents(i,:),a,b);
%         Jexpression = 0;
        
        integrandXScaled = integrandX*(sqrt(a^2*sin(theta(i))^2+b^2*cos(theta(i))^2));
        integrands(i) = integrandXScaled;
        
%         integrands(i) = integrandX;
        
        
        Jexpressions(i) = Jexpression;

    end
    
    figure; hold on;
    plot(theta,integrands)
    plot(theta,newFun(theta))
    plot(theta,newFun2(theta))
    
    figure; hold on;
    plot(theta,integrands - newFun2(theta))
    
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

function [integrandX,integrandY,Jexpression,checkX,checkY] = calcIntegrand2(theta,t,a,b)

dr = [a*(cos(theta)-1),b*sin(theta)];

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

function [integrandX] = calcIntegrand3(theta,t,a,b)

dr = [a*(cos(theta)-1),b*sin(theta)];

dydx = dr(1)/dr(2);

integrandX = (((t(1)^2-t(2)^2)*(dr(2)-dr(1)))/(1+dydx)) - ((4*t(1)*t(2)*dydx)/((1+dydx^2)^2));

%integrandY = (dr(2)/(square.^2))*Jexpression;

end

function result = calculateIntegral2(r0,a,b)
    fun = @(t)calcIntegrandVectorSimple(t,r0,a,b);

    result = integral(fun,0,1);

end

function integrands = calcIntegrandVectorSimple(t,r0,a,b)
    integrands = zeros(1,numel(t));
    
    % the theta values then go from 0 to 2pi
    theta = (2*pi).*t;
    
    tangents(:,1) = -a*sin(theta);
    tangents(:,2) = b*cos(theta);
    
    for i=1:numel(integrands)
        integrandX = calcIntegrand4(theta(i),tangents(i,:),a,b);
        integrands(i) = integrandX;
    end

end

function [integrandX] = calcIntegrand4(theta,t,a,b)

dr = [a*(cos(theta)-1),b*sin(theta)];

square = (dr(1).^2 + dr(2).^2);

tX = t(1);
tY = t(2);
dX = dr(1);
dY = dr(2);

integrandX = ((dr(1)*((a^2*sin(theta)^2+b^2*cos(theta)^2)^(1/2)))/(square.^2))*(tX*tX*(dY*dY-dX*dX) - 4*tX*tY*dX*dY + tY*tY*(dX*dX-dY*dY));

%integrandX = (a*(a^2*sin(theta)^2+b^2*cos(theta)^2)^(1/2)*(cos(theta)-1)*((a^2*sin(theta)^2*(b^2*sin(theta)^2-a^2*(cos(theta)-1)^2))/(a^2*sin(theta)^2+b^2*cos(theta)^2) + (-b^2*cos(theta)^2*(b^2*sin(theta)^2-a^2*(cos(theta)-1)^2))/(a^2*sin(theta)^2+b^2*cos(theta)^2)+(4*a^2*b^2*cos(theta)*sin(theta)^2*(cos(theta)-1))/(a^2*sin(theta)^2+b^2*cos(theta)^2)))/(b^2*sin(theta)^2+a^2*(cos(theta)-1)^2)^2;
 

%integrandY = (dr(2)/(square.^2))*Jexpression;

end