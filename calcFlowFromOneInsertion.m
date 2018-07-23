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

function flow = calcFlowFromOneInsertion(pos,insertion,insRate,R) 

disp('calculating')
disp(insertion)

flow = [0,0];

for m = -10:10
    
    sm = [0,2*pi*m*R];
    
    vectorToInsertion = pos - (insertion + sm);
    
    distToInsertion = vectorToInsertion(1)^2+vectorToInsertion(2)^2;
    
    flow = flow + vectorToInsertion./distToInsertion;
    
end

flow = flow*insRate/(2*pi);
