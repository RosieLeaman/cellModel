function [problemFlag,problemsProteinProtein,problemsLPSLPS,problemsProteinLPS,problemsProtein,problemsLPS] = checkPolygonDistances(model,tooClose)

problemFlag = 0;
problemsProteinProtein = [];
problemsLPSLPS = [];
problemsProteinLPS = [];
problemsProtein = [];
problemsLPS = [];

% first we check for protein-protein issues
index = 1;
for proteinBlob1 = 1:numel(model.proteinVertices)
    for proteinBlob2 = proteinBlob1+1:numel(model.proteinVertices)
        % we go through each vertex in proteinBlob1 and see if it's close
        % to any vertices in proteinBlob2
        for i=1:size(model.proteinVertices{proteinBlob1},1)
            x = model.proteinVertices{proteinBlob1}(i,:);
            try
                distances = model.proteinVertices{proteinBlob2} - x;
            catch
                model.proteinVertices
                assert(1==0);
            end
            distances = distances.^2;
            distances = sum(distances,2);
            
            if sum(distances < tooClose) > 0
                % at least one of these vertices is 'too close'
                problemFlag = 1;
                problemsProteinProtein(index,:) = [proteinBlob1,proteinBlob2];
                index = index + 1;
                % here we break out of going round proteinBlob1 vertices as
                % there's no point keeping going since we've already
                % decided there's an issue
                break 
            end
        end
    end
end

% next we check for lps-lps issues
index = 1;
for lpsBlob1 = 1:numel(model.lpsVertices)
    for lpsBlob2 = lpsBlob1+1:numel(model.lpsVertices)
        % we go through each vertex in lpsBlob1 and see if it's close
        % to any vertices in lpsBlob2
        for i=1:size(model.lpsVertices{lpsBlob1},1)
            x = model.lpsVertices{lpsBlob1}(i,:);
            distances = model.lpsVertices{lpsBlob2} - x;
            distances = distances.^2;
            distances = sum(distances,2);
            
            if sum(distances < tooClose) > 0
                % at least one of these vertices is 'too close'
                problemFlag = 1;
                problemsLPSLPS(index,:) = [lpsBlob1,lpsBlob2];
                index = index + 1;
                % here we break out of going round proteinBlob1 vertices as
                % there's no point keeping going since we've already
                % decided there's an issue
                break 
            end
        end
    end
end

% then we want to check every protein blob to see if its close to an LPS blob
index = 1;
for proteinBlob=1:numel(model.proteinVertices)
    for lpsBlob=1:numel(model.lpsVertices)
        % go through each protein vertex
        for i=1:size(model.proteinVertices{proteinBlob},1)
            x = model.proteinVertices{proteinBlob}(i,:);
            % find the distance from this vertex to all the LPS vertices in
            % lpsBlob
            distances = model.lpsVertices{lpsBlob} - x;
            distances = distances.^2;
            distances = sqrt(sum(distances,2)); % want to sum across rows
            
            if sum(distances < tooClose) > 0
                % at least one of these vertices is 'too close'
                problemFlag = 1;
                problemsProteinLPS(index,:) = [proteinBlob,lpsBlob];
                index = index + 1;
                % here we break out of going round proteinBlob1 vertices as
                % there's no point keeping going since we've already
                % decided there's an issue
                break 
            end
        end        
    end
end

% we then want to check for being close to self i.e. ready to split
% this is actually more complicated than I thought as neighbouring vertices
% are obviously going to be close but don't really count!
% maybe we should check if the distance round the polygon is ALSO > 1
index = 1;
finished = 0;
for proteinBlob = 1:numel(model.proteinVertices)
    for i=1:size(model.proteinVertices{proteinBlob},1)
        x = model.proteinVertices{proteinBlob}(i,:);
        for j=i+1:size(model.proteinVertices{proteinBlob},1)
            y = model.proteinVertices{proteinBlob}(j,:);
            % check if they are near as the crow flies
            if findDist(x,y) < 1
                % but also far around the polygon
                if findDistRoundPolygon(i,j,model.proteinVertices{proteinBlob}) > 2*pi
                    % but also, we need to check whether the part of the
                    % polygon between the two points in inside or outside
                    % the polygon. Outside is no problem
                    
                    % had to change this to checking a whole line, not just
                    % one point as there were issues with overhangs
                    inside = checkLineInPolygon(x,y,model.proteinVertices{proteinBlob});
                    if inside == 1
                        %disp('found a sketch')
                        problemFlag = 1;
                        problemsProtein(index) = proteinBlob;
                        index = index + 1;
                        finished = 1;
                        break
                    end
                end
            end
        end
        if finished == 1
            break
        end
    end
end

index = 1;
finished = 0;
for lpsBlob = 1:numel(model.lpsVertices)
    for i=1:size(model.lpsVertices{lpsBlob},1)
        x = model.lpsVertices{lpsBlob}(i,:);
        for j=i+1:size(model.lpsVertices{lpsBlob},1)
            y = model.lpsVertices{lpsBlob}(j,:);
            % check if they are near as the crow flies
            if findDist(x,y) < 1
                % but also far around the polygon
                if findDistRoundPolygon(i,j,model.lpsVertices{lpsBlob}) > 2
                    % but also, we need to check whether the part of the
                    % polygon between the two points in inside or outside
                    % the polygon. Outside is no problem
                    inside = checkLineInPolygon(x,y,model.lpsVertices{lpsBlob});
                    if inside == 1
                        problemFlag = 1;
                        problemsLPS(index) = lpsBlob;
                        index = index + 1;
                        finished = 1; % need an extra count so can break out of the other for loop too
                        break
                    end
                end
            end
        end
        if finished == 1
            break
        end
    end
end