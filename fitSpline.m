function [ppvalX,ppvalY] = fitSpline(x,parametrisation)
numVerticesTripled = size(x,1);
numVertices = numVerticesTripled/3;

% % first get the distances apart from each other
t = zeros(numVerticesTripled,1);
for i=2:numel(t)
    t(i) = t(i-1) + findDist(x(i,:),x(i-1,:));
end

distAround = t(numVertices+1);

assert((t(2*numVertices+1)-2*distAround)< 10e-10,'fitSpline: got different results on each loop round')

% normalise so that the whole circle is between 0 and 1 (with 0 and 1 being
% same point)
% stretch it so it's actually -1 to 2-1/numVertices
t = (t./distAround)-1;

% normalise to -1 to 2-1/numVertices
% Note this is the original numVertices, i.e. the number of vertices here
% divided by three as we copied points 3 times
% this allows us to have the middle copy of points, where the tangents
% should be good all the way round as far from end points, to be the
% section of the spline where the parameter is between 0 and 1

if parametrisation == 1
    t = linspace(-1,2-1/numVertices,numVerticesTripled);
end

% get the x and y splines separately
ppvalX = spline(t,x(:,1));
ppvalY = spline(t,x(:,2));
