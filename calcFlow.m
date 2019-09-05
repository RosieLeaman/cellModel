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

function flows = calcFlow(allPos,insertionLocsProtein,insRateProtein,insertionLocsLPS,insRateLPS,smVecsY)
%%%% new %%%%

fakeInsertions = zeros(size(insertionLocsProtein,1)*21,2);

for i=1:size(insertionLocsProtein,1)
    fakeInsertions(21*(i-1)+1:21*i,:) = insertionLocsProtein(i,:) + smVecsY;
end

% calculate the proteinFlow
%proteinFlow = zeros(size(allPos));
flows = zeros(size(allPos));

% go through each insertion point adding flow from that insertion point to
% each points flow

for i = 1:size(fakeInsertions,1)
    % distance from the point to that insertion point
    posToInsertions = allPos - fakeInsertions(i,:);
    % this is weirdly faster than norm() for very large numbers of
    % iterations
    distToInsertions = posToInsertions(:,1).^2 + posToInsertions(:,2).^2;
    
    %proteinFlow = proteinFlow + posToInsertions./distToInsertions;
    flows = flows + insRateProtein*posToInsertions./distToInsertions;
end

%proteinFlow = proteinFlow*insRateProtein;


fakeInsertions = zeros(size(insertionLocsLPS,1)*21,2);

for i=1:size(insertionLocsLPS,1)
    fakeInsertions(21*(i-1)+1:21*i,:) = insertionLocsLPS(i,:) + smVecsY;
end

% calculate the proteinFlow
%LPSFlow = zeros(size(allPos));

% go through each insertion point adding flow from that insertion point to
% each points flow

for i = 1:size(fakeInsertions,1)
    % distance from the point to that insertion point
    posToInsertions = allPos - fakeInsertions(i,:);
    % this is weirdly faster than norm() for very large numbers of
    % iterations
    distToInsertions = posToInsertions(:,1).^2 + posToInsertions(:,2).^2;
    
    %LPSFlow = LPSFlow + posToInsertions./distToInsertions;
    flows = flows + insRateLPS*posToInsertions./distToInsertions;
end

%LPSFlow = LPSFlow*insRateLPS/(2*pi);

flows = flows./(2*pi);

%flows = proteinFlow + LPSFlow;

%%%% old %%%%
% flows = zeros(size(allPos));
% 
% for vertex = 1:size(allPos,1)
%     pos = allPos(vertex,:);
% 
%     proteinFlow = [0,0];
% 
%     % calculate flow due to proteins
% 
%     for i=1:size(insertionLocsProtein,1)
%         posToInsertionVec = pos-insertionLocsProtein(i,:);
% 
%         posToInsertionVecsY = posToInsertionVec(2) + smVecsY;
% 
%         distToInsertions = posToInsertionVec(1).^2+posToInsertionVecsY.^2;
% 
%         flow = [sum(posToInsertionVec(1)./distToInsertions),sum(posToInsertionVecsY./distToInsertions)];
% 
%         proteinFlow = proteinFlow + flow;
%     end
% 
%     proteinFlow = proteinFlow*insRateProtein/(2*pi);
% 
% %     % calculate flow due to LPS
% % 
%      LPSflow = [0,0];
% % 
% %     for i=1:size(insertionLocsLPS,1)
% %         posToInsertionVec = pos-insertionLocsLPS(i,:);
% % 
% %         posToInsertionVecsY = posToInsertionVec(2) + smVecsY;
% % 
% %         distToInsertions = posToInsertionVec(1).^2+posToInsertionVecsY.^2;
% % 
% %         flow = [sum(posToInsertionVec(1)./distToInsertions),sum(posToInsertionVecsY./distToInsertions)];
% % 
% %         LPSflow = LPSflow + flow;
% %     end
% % 
%     LPSflow = LPSflow*insRateLPS/(2*pi);
% 
%     flows(vertex,:) = proteinFlow + LPSflow;
%     flows(vertex,:) = proteinFlow;
% end