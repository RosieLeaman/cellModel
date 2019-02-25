% finds the tangent at x1

function [normal,tangent] = findTangentQuadratic(x0,x1,x2,unit)

% fit a quadratic to the three points
p = polyfit([x0(1),x1(1),x2(1)],[x0(2),x1(2),x2(2)],2);

% work out the gradient of the quadratic at the point x1

m = 2*p(1)*x1(1) + p(2);

% work out what c should be

c = x1(2) - x1(1)*m;

% work out what the value should be at x2(1)

y2 = x2(1)*m + c;

% then work out a vector by subtracting the two

tangent = [x2(1),y2] - x1;

% check that the tangent points clockwise the whole way round
%

if x1(1) >= 0 && tangent(1) < 0 || x1(1) < 0 && tangent(1) > 0
    % if we are in top half plane then the x tangent component should be
    % pointing right, i.e. positive. Similarly if we are in the bottom half
    % plane the x tangent should be pointing left, i.e. negative.
    % if this is not true we should make the tangent negative
    
    tangent = -tangent;
end

% make it a unit vector if desired
if unit == 1
    tangentNorm = sqrt(dot(tangent,tangent));
    tangent = tangent./tangentNorm;
end

% construct the outward pointing normal
normal = [-tangent(2),tangent(1)];

