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

function flow = calcFlow(pos,insertionLocsProtein,insRateProtein,insertionLocsLPS,insRateLPS,membraneCircumference,flag)

smVecsY = ((-10:10)')*membraneCircumference;

proteinFlow = [0,0];

% calculate flow due to proteins

for i=1:size(insertionLocsProtein,1)
    %if sum((pos-insertionLocsProtein(i,:)).^2) < 22500
    %nextFlow = calcFlowFromOneInsertion(pos,insertionLocsProtein(i,:),membraneCircumference);

    posToInsertionVec = pos-insertionLocsProtein(i,:);

    posToInsertionVecsY = posToInsertionVec(2) + smVecsY;

    distToInsertions = posToInsertionVec(1).^2+posToInsertionVecsY.^2;

    flow = [sum(posToInsertionVec(1)./distToInsertions),sum(posToInsertionVecsY./distToInsertions)];

    proteinFlow = proteinFlow + flow;
    %end
end

proteinFlow = proteinFlow*insRateProtein/(2*pi);

% calculate flow due to LPS

LPSflow = [0,0];

for i=1:size(insertionLocsLPS,1)
    nextFlow = calcFlowFromOneInsertion(pos,insertionLocsLPS(i,:),membraneCircumference);
%     if flag == 1
%         disp(['flow due to lptd ',num2str(i)])
%         nextFlow
%         pos
%         insertionLocsLPS
%         insRateLPS
%     end
    LPSflow = LPSflow + nextFlow;
end

LPSflow = LPSflow*insRateLPS/(2*pi);

flow = proteinFlow + LPSflow;

end

function flow = calcFlowFromOneInsertion(pos,insertion,membraneCircumference) 

smVecsY = ((-10:10)')*membraneCircumference;

posToInsertionVec = pos-insertion;

posToInsertionVecsY = posToInsertionVec(2) + smVecsY;

distToInsertions = posToInsertionVec(1).^2+posToInsertionVecsY.^2;

flow = [sum(posToInsertionVec(1)./distToInsertions),sum(posToInsertionVecsY./distToInsertions)];

end