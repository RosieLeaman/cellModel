% this function tests the effectiveness of the spline at estimating the
% tangent

function [errorX,errorY] = testSplineTangent(numVertices,a,b,plotYes)

format long
% if there are no vertices input we just use a default value
if nargin < 1
    numVertices = 1000;
end

if nargin <= 1
    % no ellipse parameters specified, use defaults
    a = 4; % x axis stretch
    b = 2; % y axis stretch
end

[points,angles] = findVerticesNewMaterialEllipse([0,0],numVertices,a,b);
%[points,angles] = findVerticesNewMaterialEllipseWithError([0,0],numVertices,a,b,10e-4);

 % to get nice estimates of the tangent at the 'beginning'/'end' of the
% polygon list, we extend the points list by tripling it. Then the
% central 'points' list has good estimates of the tangent
pointsTripled = [points;points;points];

% we fit a spline to these points
[splineX,splineY] = fitSpline(pointsTripled);

% from the spline we calculate a tangent and normals at the polygon
% points. We actually only need the normals which are used to move the
% points along the normal direction
[wholeTangents,wholeNormals] = findTangentFromSplines(splineX.breaks,splineX,splineY,1);

% the above gives us normals and tangents for each point going around
% the polygon three times. Just select out the middle bit that we
% actually want here
normals = wholeNormals(numVertices+1:numVertices*2,:);
tangents = wholeTangents(numVertices+1:numVertices*2,:);

% test the tangent is close to accurate
% the analytical unit tangent should be for an ellipse
% t = (-(a/b)*y,(b/a)*x). With unit tangent being t./||t||
x = linspace(0,2*pi*(1-1/numVertices),numVertices);
bigX = linspace(0,2*pi*(1-1/(2*numVertices)),2*numVertices);

testT = zeros(2*numVertices,1);
testIndex = 1;
for i=1:numVertices
    testT(testIndex) = splineX.breaks(numVertices+i);
    testT(testIndex+1) = (splineX.breaks(numVertices+i)+splineX.breaks(numVertices+i+1))/2;
    testIndex = testIndex + 2;
end

testX = ppval(splineX,testT);
testY = ppval(splineY,testT);

thetas = atan2(a*testY,b*testX)';
thetasOriginal = atan2(a*points(:,2),b*points(:,1))';

newTangents = findTangentFromSplines(testT',splineX,splineY,1);

errorVecX = abs(-a*sin(thetas)./sqrt(a*a*sin(thetas).^2+b*b*cos(thetas).^2)-newTangents(:,1)');
errorVecY = abs(b*cos(thetas)./sqrt(a*a*sin(thetas).^2+b*b*cos(thetas).^2)-newTangents(:,2)');

errorVecOriginalX = abs(-a*sin(thetasOriginal)./sqrt(a*a*sin(thetasOriginal).^2+b*b*cos(thetasOriginal).^2)-tangents(:,1)');
errorVecOriginalY = abs(b*cos(thetasOriginal)./sqrt(a*a*sin(thetasOriginal).^2+b*b*cos(thetasOriginal).^2)-tangents(:,2)');

errorX = mean(errorVecX);
errorY = mean(errorVecY);

% errorX = mean(errorVecOriginalX);
% errorY = mean(errorVecOriginalY);


% figure;
% subplot(1,2,1); hold on;
% plot(abs(-a*sin(thetasOriginal)./sqrt(a*a*sin(thetasOriginal).^2+b*b*cos(thetasOriginal).^2)-tangents(:,1)'),'x-')
% title('x')
% 
% subplot(1,2,2); hold on;
% plot(abs(b*cos(thetasOriginal)./sqrt(a*a*sin(thetasOriginal).^2+b*b*cos(thetasOriginal).^2)-tangents(:,2)'),'x-')
% title('y')

if plotYes == 1
    figure;
    subplot(1,2,1); hold on;
    plot(errorVecX,'-')
    title('x')

    subplot(1,2,2); hold on;
    plot(errorVecY,'-')
    title('y')
end