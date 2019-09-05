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

numBAM = size(insertionLocsProtein,1);

fakeInsertions = zeros(numBAM*21,2);

for i=1:size(insertionLocsProtein,1)
    fakeInsertions(21*(i-1)+1:21*i,:) = insertionLocsProtein(i,:) + smVecsY;
end

% calculate the proteinFlow
flows = zeros(size(allPos));

% go through each insertion point adding flow from that insertion point to
% each points flow

for i = 1:size(fakeInsertions,1)
    % distance from the point to that insertion point
    posToInsertions = allPos - fakeInsertions(i,:);
    % this is weirdly faster than norm() for very large numbers of
    % iterations
    distToInsertions = (posToInsertions(:,1).^2 + posToInsertions(:,2).^2);
    
    invDists = 1./distToInsertions;
    
    newFlow = insRateProtein*posToInsertions.*invDists;
    
    flows = flows + newFlow;
end

% do same for LPS
fakeInsertions = zeros(size(insertionLocsLPS,1)*21,2);

for i=1:size(insertionLocsLPS,1)
    fakeInsertions(21*(i-1)+1:21*i,:) = insertionLocsLPS(i,:) + smVecsY;
end

% go through each insertion point adding flow from that insertion point to
% each points flow

for i = 1:size(fakeInsertions,1)
    % distance from the point to that insertion point
    posToInsertions = allPos - fakeInsertions(i,:);
    % this is weirdly faster than norm() for very large numbers of
    % iterations
    distToInsertions = posToInsertions(:,1).^2 + posToInsertions(:,2).^2;
    
    invDists = 1./distToInsertions;
    
    newFlow = insRateLPS*posToInsertions.*invDists;
    
    flows = flows + newFlow;
end

flows = flows./(2*pi);



% numBAM = size(insertionLocsProtein,1);
% 
% fakeInsertions = zeros(numBAM*21,2);
% 
% for i=1:size(insertionLocsProtein,1)
%     fakeInsertions(21*(i-1)+1:21*i,:) = insertionLocsProtein(i,:) + smVecsY;
% end
% 
% % calculate the proteinFlow
% flows = zeros(size(allPos));
% 
% % go through each insertion point adding flow from that insertion point to
% % each points flow
% 
% for i = 1:size(fakeInsertions,1)
%     % distance from the point to that insertion point
%     posToInsertions = allPos - fakeInsertions(i,:);
%     % this is weirdly faster than norm() for very large numbers of
%     % iterations
%     distToInsertions = posToInsertions(:,1).^2 + posToInsertions(:,2).^2;
%     
%     flows = flows + insRateProtein*posToInsertions./distToInsertions;
% end
% 
% % do same for LPS
% fakeInsertions = zeros(size(insertionLocsLPS,1)*21,2);
% 
% for i=1:size(insertionLocsLPS,1)
%     fakeInsertions(21*(i-1)+1:21*i,:) = insertionLocsLPS(i,:) + smVecsY;
% end
% 
% % go through each insertion point adding flow from that insertion point to
% % each points flow
% 
% for i = 1:size(fakeInsertions,1)
%     % distance from the point to that insertion point
%     posToInsertions = allPos - fakeInsertions(i,:);
%     % this is weirdly faster than norm() for very large numbers of
%     % iterations
%     distToInsertions = posToInsertions(:,1).^2 + posToInsertions(:,2).^2;
%     
%     flows = flows + insRateLPS*posToInsertions./distToInsertions;
% end
% 
% flows = flows./(2*pi);
% 
