function [tangents,normals] = findTangentFromSplines(t,splineX,splineY,unit)

% we now need to 'differentiate' these to get dx/dt and dy/dt

splineXdiff = fnder(splineX);
splineYdiff = fnder(splineY);

xdiffs = ppval(splineXdiff,t);
ydiffs = ppval(splineYdiff,t);

tangents = [xdiffs',ydiffs'];

if unit == 1
    for i=1:size(tangents,1)
        tangents(i,:) = tangents(i,:)./norm(tangents(i,:));
    end
end

normals = zeros(size(tangents));
for i=1:size(tangents,1)
    % this is the OUTWARD POINTING NORMAL
    normals(i,:) = [tangents(i,2),-tangents(i,1)];
end