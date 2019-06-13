function [normals,tangents,ppvalX,ppvalY] = findTangentSpline(x,unit)

% % first get the distances apart from each other
% t = zeros(size(x,1),1);
% for i=2:numel(t)
%     t(i) = t(i-1) + findDist(x(i,:),x(i-1,:));
% end
% 
% % DO NORMALISE. THIS IS EASIER.
% % normalise to [0,1]
% t = t./(max(t));

% % do it with theta
% t = zeros(size(x,1),1);
% for i=2:numel(t)
%     t(i) = atan2(x(i,2),x(i,1));
%     if t(i) <= 0
%         t(i)= t(i) + 2*pi;
%     end
% end
% t = t./(max(t));

% % do it with theta
% t = zeros(size(x,1),1);
% a = 0;
% for i=2:(numel(t))
%     t(i) = atan2(x(i,2),x(i,1));
%     if t(i) < 0
%         t(i)= t(i) + 2*pi;
%     end
%     t(i) = t(i) + a*2*pi;
%     if i > numel(t)/3
%         t(i) = t(i) + 2*pi;
%     end
%     if i > 2*numel(t)/3
%         t(i) = t(i) + 2*pi;
%     end
% end
% t = (t-2*pi)./(2*pi);
% 
% 
% maxVal = 6*pi-(6*pi/numel(t));
% t = linspace(0,maxVal,numel(t));
% t = (t - 2*pi)./(2*pi);
% 
% t = zeros(size(x,1),1);
% for i=2:numel(t)
%     t(i) = t(i-1) + findDist(x(i,:),x(i-1,:));
% end
% 
% % DO NORMALISE. THIS IS EASIER.
% % normalise to [0,1]
% t = (t./(max(t)))*3-1;
% 
% figure;plot(linspace(0,1,numel(t)),t,'x')

% normalise to -1 to 2-1/numVertices
% Note this is the original numVertices, i.e. the number of vertices here
% divided by three as we copied points 3 times
% this allows us to have the middle copy of points, where the tangents
% should be good all the way round as far from end points, to be the
% section of the spline where the parameter is between 0 and 1

numVerticesTripled = size(x,1);
numVertices = numVerticesTripled/3;
t = linspace(-1,2-1/numVertices,numVerticesTripled);

% get the x and y splines separately
ppvalX = spline(t,x(:,1));
ppvalY = spline(t,x(:,2));

ppx = ppval(ppvalX,t);
ppy = ppval(ppvalY,t);

%figure;plot(ppx,ppy,'kx-')

% we now need to 'differentiate' these to get dx/dt and dy/dt

ppvalXdiff = fnder(ppvalX);
ppvalYdiff = fnder(ppvalY);

xdiffs = ppval(ppvalXdiff,t);
ydiffs = ppval(ppvalYdiff,t);

tangents = [xdiffs',ydiffs'];

if unit == 1
    for i=1:size(tangents,1)
        tangents(i,:) = tangents(i,:)./norm(tangents(i,:));
    end
end

normals = zeros(size(x,1),2);

for i=1:size(tangents,1)
    % this is the OUTWARD POINTING NORMAL
    normals(i,:) = [-tangents(i,2),tangents(i,1)];
end