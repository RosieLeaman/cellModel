% checks to see if any polygons are `too close' where tooClose is a number
% defining how close is too close, and if they are joins them into one
% polygon. Issues are resolved once per iteration.

function [problems,newModel] = checkPolygonDistances2(model,tooClose,plotYes)

% make a new model
newModel = model;
problems = 0;

% first we check for whether two different protein regions are close
% together, and should be joined into one.
foundProblem = 0;

for proteinBlob1 = 1:numel(model.proteinVertices)
    for proteinBlob2 = proteinBlob1+1:numel(model.proteinVertices)
        % we go through each vertex in proteinBlob1 and see if it's close
        % to any vertices in proteinBlob2
        for i=1:size(model.proteinVertices{proteinBlob1},1)
            x = model.proteinVertices{proteinBlob1}(i,:);
            
            % subtract the point x from each vertex in proteinBlob2
            distances = model.proteinVertices{proteinBlob2} - x;
            distances = distances.^2; % squared distance
            distances = sum(distances,2); % sum the x and y squares for each point
            
            % this check makes a logical vector with a 1 if the distance in
            % that index is < tooClose. If the vector is summed and the sum
            % is > 0 then at least one point must be too close
            if sum(distances < tooClose) > 0
                foundProblem = 1;
                % escape the loop and resolve this issue
                break 
            end
        end
        if foundProblem == 1
            break
        end
    end
    if foundProblem == 1
        break
    end
end

if foundProblem == 1
    problems = 1;
    if plotYes == 1
        figure;
        visualiseSimple(model)
    end
    
    %disp(['an issue was found between polygons ',num2str(proteinBlob1),' and ',num2str(proteinBlob2)])
    
    % an issue was found, resolve it.
    newVertices = resolveTwoSame(model.proteinVertices{proteinBlob1},model.proteinVertices{proteinBlob2},plotYes);
    
    % remove one of the polygons
    newModel.proteinVertices(proteinBlob2) = [];
    
    % replace the other polygon with the new vertices
    newModel.proteinVertices{proteinBlob1} = newVertices;
end

% 
% % next we check for lps-lps issues
% index = 1;
% for lpsBlob1 = 1:numel(model.lpsVertices)
%     for lpsBlob2 = lpsBlob1+1:numel(model.lpsVertices)
%         % we go through each vertex in lpsBlob1 and see if it's close
%         % to any vertices in lpsBlob2
%         for i=1:size(model.lpsVertices{lpsBlob1},1)
%             x = model.lpsVertices{lpsBlob1}(i,:);
%             distances = model.lpsVertices{lpsBlob2} - x;
%             distances = distances.^2;
%             distances = sum(distances,2);
%             
%             if sum(distances < tooClose) > 0
%                 % at least one of these vertices is 'too close'
%                 problemFlag = 1;
%                 problemsLPSLPS(index,:) = [lpsBlob1,lpsBlob2];
%                 index = index + 1;
%                 % here we break out of going round proteinBlob1 vertices as
%                 % there's no point keeping going since we've already
%                 % decided there's an issue
%                 break 
%             end
%         end
%     end
% end
% 
% % then we want to check every protein blob to see if its close to an LPS blob
% index = 1;
% for proteinBlob=1:numel(model.proteinVertices)
%     for lpsBlob=1:numel(model.lpsVertices)
%         % go through each protein vertex
%         for i=1:size(model.proteinVertices{proteinBlob},1)
%             x = model.proteinVertices{proteinBlob}(i,:);
%             % find the distance from this vertex to all the LPS vertices in
%             % lpsBlob
%             distances = model.lpsVertices{lpsBlob} - x;
%             distances = distances.^2;
%             distances = sqrt(sum(distances,2)); % want to sum across rows
%             
%             if sum(distances < tooClose) > 0
%                 % at least one of these vertices is 'too close'
%                 problemFlag = 1;
%                 problemsProteinLPS(index,:) = [proteinBlob,lpsBlob];
%                 index = index + 1;
%                 % here we break out of going round proteinBlob1 vertices as
%                 % there's no point keeping going since we've already
%                 % decided there's an issue
%                 break 
%             end
%         end        
%     end
% end
% 
% % we then want to check for being close to self i.e. ready to split
% % this is actually more complicated than I thought as neighbouring vertices
% % are obviously going to be close but don't really count!
% % maybe we should check if the distance round the polygon is ALSO > 1
% index = 1;
% finished = 0;
% for proteinBlob = 1:numel(model.proteinVertices)
%     for i=1:size(model.proteinVertices{proteinBlob},1)
%         x = model.proteinVertices{proteinBlob}(i,:);
%         for j=i+1:size(model.proteinVertices{proteinBlob},1)
%             y = model.proteinVertices{proteinBlob}(j,:);
%             % check if they are near as the crow flies
%             if findDist(x,y) < 1
%                 % but also far around the polygon
%                 if findDistRoundPolygon(i,j,model.proteinVertices{proteinBlob}) > 2*pi
%                     % but also, we need to check whether the part of the
%                     % polygon between the two points in inside or outside
%                     % the polygon. Outside is no problem
%                     
%                     % had to change this to checking a whole line, not just
%                     % one point as there were issues with overhangs
%                     inside = checkLineInPolygon(x,y,model.proteinVertices{proteinBlob});
%                     if inside == 1
%                         %disp('found a sketch')
%                         problemFlag = 1;
%                         problemsProtein(index) = proteinBlob;
%                         index = index + 1;
%                         finished = 1;
%                         break
%                     end
%                 end
%             end
%         end
%         if finished == 1
%             break
%         end
%     end
% end
% 
% index = 1;
% finished = 0;
% for lpsBlob = 1:numel(model.lpsVertices)
%     for i=1:size(model.lpsVertices{lpsBlob},1)
%         x = model.lpsVertices{lpsBlob}(i,:);
%         for j=i+1:size(model.lpsVertices{lpsBlob},1)
%             y = model.lpsVertices{lpsBlob}(j,:);
%             % check if they are near as the crow flies
%             if findDist(x,y) < 1
%                 % but also far around the polygon
%                 if findDistRoundPolygon(i,j,model.lpsVertices{lpsBlob}) > 2
%                     % but also, we need to check whether the part of the
%                     % polygon between the two points in inside or outside
%                     % the polygon. Outside is no problem
%                     inside = checkLineInPolygon(x,y,model.lpsVertices{lpsBlob});
%                     if inside == 1
%                         problemFlag = 1;
%                         problemsLPS(index) = lpsBlob;
%                         index = index + 1;
%                         finished = 1; % need an extra count so can break out of the other for loop too
%                         break
%                     end
%                 end
%             end
%         end
%         if finished == 1
%             break
%         end
%     end
% end