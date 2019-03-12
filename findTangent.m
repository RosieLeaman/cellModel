% This function takes in three points which should be spaced along the
% curve x0 - x1 - x2 and estimates the normal and the tangent at the point
% x1. The normal returned is the OUTWARD POINTING normal
% if we have the tangent at x1, [t1,t2] then the normal is [-t2,t1]
% unit should be 1 if unit vectors are desired and zero otherwise

% OUTPUTS:
% normal; a 1x2 row vector which is the outward pointing normal to the curve at x1
% tangent; a 1x2 row vector which is the tangent to the curve at x1

function [normal,tangent] = findTangent(x1,x2,unit)

tangent = 0.5*(x2-x1);

if unit == 1
    tangentNorm = sqrt(dot(tangent,tangent));
    tangent = tangent./tangentNorm;
end

normal = [tangent(2),-tangent(1)];

