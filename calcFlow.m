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

function flows = calcFlow(allPos,insertionLocsProtein,insRateProtein,insertionLocsLPS,insRateLPS,membraneCircumference,smVecsY)

flows = zeros(size(allPos));

for vertex = 1:size(allPos,1)
    pos = allPos(vertex,:);

    proteinFlow = [0,0];

    % calculate flow due to proteins

    for i=1:size(insertionLocsProtein,1)
        posToInsertionVec = pos-insertionLocsProtein(i,:);

        posToInsertionVecsY = posToInsertionVec(2) + smVecsY;

        distToInsertions = posToInsertionVec(1).^2+posToInsertionVecsY.^2;

        flow = [sum(posToInsertionVec(1)./distToInsertions),sum(posToInsertionVecsY./distToInsertions)];

        proteinFlow = proteinFlow + flow;
    end

    proteinFlow = proteinFlow*insRateProtein/(2*pi);

    % calculate flow due to LPS

    LPSflow = [0,0];

    for i=1:size(insertionLocsLPS,1)
        posToInsertionVec = pos-insertionLocsLPS(i,:);

        posToInsertionVecsY = posToInsertionVec(2) + smVecsY;

        distToInsertions = posToInsertionVec(1).^2+posToInsertionVecsY.^2;

        flow = [sum(posToInsertionVec(1)./distToInsertions),sum(posToInsertionVecsY./distToInsertions)];

        LPSflow = LPSflow + flow;
    end

    LPSflow = LPSflow*insRateLPS/(2*pi);

    flows(vertex,:) = proteinFlow + LPSflow;
end