% finds the tangent at x1

function [normal,tangent] = findTangentQuadratic(x0,x1,x2,unit)

% check whether our points are really vertically aligned (there are repeated x's)
%If they are then the quadratic fitting y to x will be really badly conditioned
% and won't work so it would be better to fit x to y

verticalFlag = 0;

warning('') % Clear last warning message
    

if abs(x0(1) - x1(1)) < 10e-5 || abs(x1(1) - x2(1)) < 10e-5 || abs(x0(1) - x2(1)) < 10e-5
    % in this case we want to fit x to y not y to x
    verticalFlag = 1;
    x0 = fliplr(x0);
    x1 = fliplr(x1);
    x2 = fliplr(x2);
end

% fit a quadratic to the three points
p = polyfit([x0(1),x1(1),x2(1)],[x0(2),x1(2),x2(2)],2);

[warnMsg, warnId] = lastwarn;
if ~isempty(warnMsg)
    disp(verticalFlag)
    disp(x0)
    disp(x1)
    disp(x2)
    disp(p)
    
    abs(x0(1) - x2(1))
    abs(x0(1) - x2(1)) < 10e-12
    error('all wrong')
end

% work out the gradient of the quadratic at the point x1

m = 2*p(1)*x1(1) + p(2);

% work out what c should be

c = x1(2) - x1(1)*m;

% work out what the value should be at x2(1)

y2 = x2(1)*m + c;

% then work out a vector by subtracting the two

tangent = [x2(1),y2] - x1;

if verticalFlag == 1
    % we want to swap back the tangent
    tangent = fliplr(tangent);
end

% check that the tangent points clockwise the whole way round
% has to be <= 0 for the tangent and >= 0 for the x to make sure we catch
% when y=0 or x=0 cases

% removed this bit for now as can't quite work out how to make it work and
% doesn't appear to be necessary as we are going anti clockwise around the
% circle anyway in a cyclical way so the tangent points the right way
% anyway

% if x1(2) >= 0 && tangent(1) < 0 || x1(2) < 0 && tangent(1) >= 0
%     % if we are in top half plane then the x tangent component should be
%     % pointing right, i.e. positive. Similarly if we are in the bottom half
%     % plane the x tangent should be pointing left, i.e. negative.
%     % if this is not true we should make the tangent negative
%     
%     tangent = -tangent;
% end

% make it a unit vector if desired
if unit == 1
    tangentNorm = norm(tangent);
    tangent = tangent./tangentNorm;
end

% construct the outward pointing normal
normal = [tangent(2),-tangent(1)];

