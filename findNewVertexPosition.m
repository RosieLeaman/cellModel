function newVertex = findNewVertexPosition(vertex,flow,dt,membraneCircumference)

newVertex = vertex + flow*dt;

% % check to see if things have fallen off the edge
% if newVertex(2) > membraneCircumference % top edge
%     newVertex(2) = newVertex(2) - membraneCircumference;
% elseif newVertex(2) < 0 % bottom edge
%     newVertex(2) = newVertex(2) + membraneCircumference;  
% end
