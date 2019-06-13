function [ppvalX,ppvalY] = fitSpline(x)

% % first get the distances apart from each other
% t = zeros(size(x,1),1);
% for i=2:numel(t)
%     t(i) = t(i-1) + findDist(x(i,:),x(i-1,:));
% end
% 
% % DO NORMALISE. THIS IS EASIER.
% % normalise to [0,1]
% t = t./(max(t));

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
