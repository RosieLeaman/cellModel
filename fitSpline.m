% NOTE: this fitSpline function is designed to work with surface tension,
% so it expects a list of points that define a polygon, where the points
% have been tripled three times so they go like
% [[x0 y0];
%    ...;
%  [xN yN];
%  [x0 y0];
%    ...;
%  [xN yN];
%  [x0 y0];
%    ...;
%  [xN yN]];

% it fits an x and y spline to these points and returns the spline object.
% The x and y splines are parametrised by the distance between points

% INPUTS
% x; an 3mxn matrix of row vectors as defined above

% OUTPUTS
% ppvalX; a spline object which passes through the x values provided
% ppvalY; a spline object which passes through the y values provided
function [ppvalX,ppvalY] = fitSpline(x)

% we find out how many vertices there were originally
numVerticesTripled = size(x,1);
numVertices = numVerticesTripled/3; 
% numVertices can be expected to be a whole number as x should be three
% copies of points

% % first get the distances apart from each other
t = zeros(numVerticesTripled,1);
for i=2:numel(t)
    t(i) = t(i-1) + findDist(x(i,:),x(i-1,:));
end

distAround = t(numVertices+1);

% check that we got the same distance around each time
assert((t(2*numVertices+1)-2*distAround)< 10e-10,'fitSpline: got different results on each loop round')

% normalise so that the whole circle is between 0 and 1 (with 0 and 1 being
% same point)
% stretch it so it's actually -1 to 2-1/numVertices
t = (t./distAround)-1;

% get the x and y splines separately
ppvalX = spline(t,x(:,1));
ppvalY = spline(t,x(:,2));
