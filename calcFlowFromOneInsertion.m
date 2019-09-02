% this function will calculate the flow field at a point pos caused by the
% flow from an insertion point at point insertion. The formula is given in
% Ursell et. al. in the section Computational modeling of OM 
% insertion recapitulates puncta behavior (equation for j(x,x_i))

% INPUTS:
% pos; 1x2 row vector, position where we want to calculate the flow
% insertion; 1x2 row vector, position where the insertion point is
% R; scalar, radius of cylinder
% insRate; insertion rate of the insertion point

% OUTPUTS:
% flow; 1x2 row vector, flow felt in the x and y directions at the point
% pos

function flow = calcFlowFromOneInsertion(pos,insertion,insRate,membraneCircumference) 

%flow = [0,0];

smVecs = [zeros(21,1),(-10:10)']*membraneCircumference;
posVecs = repmat(pos,21,1); % repeat the vector here
insertionVecs = repmat(insertion,21,1);

vectorsToInsertions = posVecs - insertionVecs + smVecs;

distToInsertions = sum((vectorsToInsertions.^2),2);

flow = sum((vectorsToInsertions./distToInsertions),1)*insRate/(2*pi);
