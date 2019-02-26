function [normals,tangents,ppvalX,ppvalY] = findTangentSpline(x,unit)

% first get the distances apart from each other
t = zeros(size(x,1),1);
for i=2:numel(t)
    t(i) = t(i-1) + findDist(x(i,:),x(i-1,:));
end

% normalise to [0,1]
t = t./(max(t));

% get the x and y splines separately
ppvalX = spline(t,x(:,1));
ppvalY = spline(t,x(:,2));

ppx = ppval(ppvalX,t);
ppy = ppval(ppvalY,t);

figure;plot(ppx,ppy,'kx-')

% we now need to 'differentiate' these to get dx/dt and dy/dt

ppvalXdiff = fnder(ppvalX);
ppvalYdiff = fnder(ppvalY);

xdiffs = ppval(ppvalXdiff,t);
ydiffs = ppval(ppvalYdiff,t);

tangents = [xdiffs,ydiffs];

if unit == 1
    for i=1:size(tangents,1)
        tangents(i,:) = tangents(i,:)/norm(tangents(i,:));
    end
end

normals = zeros(size(x,1),2);

for i=1:size(tangents,1)
    normals(i,:) = [tangents(i,2),-tangents(i,1)];
end