% this function calculates the total flow felt at a point pos from all the
% insertion points across the cylinder
% it relies on calcFlowFromOneInsertion to run
% the formula being calculated is given in Ursell et. al. 2012 in the
% section Computational modeling of OM 
% insertion recapitulates puncta behavior (equation for J(x,x_i))

% INPUTS
% pos; 1x2 row vector where we want to calculate the flow
% insertionLocsProtein; nx2 matrix where each row is the x,y co-ordinates
% of an insertion point of protein
% insRateProtein; scalar, rate of insertion of protein
% insertionLocsLPS; mx2 matrix where each row is the x-y co-ordinates of an
% insertion point of LPS
% insRateLPS; scalar, rate of insertion of LPS
% R; scalar, radius of cylinder

% OUTPUTS
% flow; 1x2 row vector which is the flow felt at pos caused by the
% insertions

function flow = calcFlow(pos,insertionLocsProtein,insRateProtein,insertionLocsLPS,insRateLPS,R)

flow = [0,0];

% calculate flow due to proteins

for i=1:size(insertionLocsProtein,1)
    disp(['calc for protein',num2str(i)])
    insertionLocsProtein(i)
    flow = flow + calcFlowFromOneInsertion(pos,insertionLocsProtein(i,:),insRateProtein,R);
end

% calculate flow due to LPS

for i=1:size(insertionLocsLPS,1)
    disp('calc for LPS')
    flow = flow + calcFlowFromOneInsertion(pos,insertionLocsLPS(i,:),insRateLPS,R);
end