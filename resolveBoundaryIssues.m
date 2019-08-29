% sort out issues with boundaries being too close
% This returns a new model which is the same but with changed

% This does not currently take in plotYes, it just doesn't plot (the zeros
% in the calls to resolvexx or splitBlob need to be replaced by 1/plotYes)

function newModel = resolveBoundaryIssues(model,problemsProteinProtein,problemsLPSLPS,problemsProteinLPS,problemsProtein,problemsLPS)

% make a copy of model

newModel = model;

% identify problems

% this is done before passing to here

toBeChangedProtein = [];
indexProtein = 1;

for i=1:size(problemsProteinProtein,1)

    % this is a two the same so one of the polygons will be deleted later

    % find the new vertices
    j = problemsProteinProtein(i,1);
    k = problemsProteinProtein(i,2);
    newVertices = resolveTwoSame(newModel.proteinVertices{j},newModel.proteinVertices{k},0);
    
    % change the existing vertices
    
    newModel.proteinVertices{j} = newVertices;
    newModel.proteinVertices{k} = newVertices;
    
    if numel(toBeChangedProtein) == 0
        toBeChangedProtein(indexProtein) = j;
        toBeChangedProtein(indexProtein+1) = k;
        indexProtein = indexProtein + 2;
    else
        if ~ismember(toBeChangedProtein,j)
            toBeChangedProtein(indexProtein) = j;
            indexProtein = indexProtein + 1;
        end
        if ~ismember(toBeChangedProtein,k)
            toBeChangedProtein(indexProtein+1) = k;
            indexProtein = indexProtein + 1;
        end       
    end
    
     
end

% lps-lps issues

toBeChangedLPS = [];
indexLPS = 1;

for i=1:size(problemsLPSLPS,1)

    % this is a two the same so one of the polygons will be deleted later

    % find the new vertices
    j = problemsLPSLPS(i,1);
    k = problemsLPSLPS(i,2);
    newVertices = resolveTwoSame(newModel.lpsVertices{j},newModel.lpsVertices{k},0);
    
    % change the existing vertices
    
    newModel.lpsVertices{j} = newVertices;
    newModel.lpsVertices{k} = newVertices;

    if numel(toBeChangedLPS) == 0
        toBeChangedLPS(indexLPS) = j;
        toBeChangedLPS(indexLPS+1) = k;
        indexLPS = indexLPS + 2;
    else
        if ~ismember(toBeChangedLPS,j)
            toBeChangedLPS(indexLPS) = j;
            indexLPS = indexLPS + 1;
        end
        if ~ismember(toBeChangedLPS,k)
            toBeChangedLPS(indexLPS+1) = k;
            indexLPS = indexLPS + 1;
        end
    end

end

% lps-protein issues

haveToBeDeletedProtein = [];
indexProtein = 1;

haveToBeDeletedLPS = [];
indexLPS = 1;

for i=1:size(problemsProteinLPS,1)

    % this is a two the same so one of the polygons will be deleted later

    % find the new vertices
    j = problemsProteinLPS(i,1);
    k = problemsProteinLPS(i,2);
    [newVertices,insideIndex] = resolveTwoDifferent(newModel.proteinVertices{j},newModel.lpsVertices{k},0);
    
    % change the existing vertices
    
    newModel.proteinVertices{j} = newVertices;
    newModel.lpsVertices{k} = newVertices;

    if insideIndex == 1
        haveToBeDeletedProtein(indexProtein) = j;
        indexProtein = indexProtein + 1;
    else
        haveToBeDeletedLPS(indexLPS) = k;
        indexLPS = indexLPS + 1;
    end

end

% singles don't have to be changed, but new vertices have to be added.
newVerticesProtein = {};

index = 1;

for i=1:numel(problemsProtein)
    j = problemsProtein(i);
    
    [newVertices1,newVertices2] = splitBlob(newModel.proteinVertices{j},0);
    
    % set the new vertices to replace the old ones
    
    newModel.proteinVertices{j} = newVertices1;
    
    % save the other new vertices to add in a minute (can't add now or
    % it'll wreck our indexing)
    
    newVerticesProtein{index} = newVertices2;
    index = index + 1;
end

% do the same for LPS

newVerticesLPS = {};

index = 1;
for i=1:numel(problemsLPS)
    j = problemsLPS(i);
    
    [newVertices1,newVertices2] = splitBlob(newModel.lpsVertices{j},0);
    
    % set the new vertices to replace the old ones
    
    newModel.lpsVertices{j} = newVertices1;
    
    % save the other new vertices to add in a minute (can't add now or
    % it'll wreck our indexing)
    
    newVerticesLPS{index} = newVertices2;
    index = index + 1;
end

% resolve all those that need deleting and adding. This has to be done IN A
% PARTICULAR ORDER

% first we go through and check for copies and add the indices of copies to
% haveToBeDeleted

% first remove things from consideration that are going to be deleted
% anyway

toBeChangedProtein = setdiff(toBeChangedProtein,haveToBeDeletedProtein);

index = numel(haveToBeDeletedProtein)+1;
for i=1:numel(toBeChangedProtein)
    vertices1 = newModel.proteinVertices{toBeChangedProtein(i)};
    for j=i+1:numel(toBeChangedProtein)
        vertices2 = newModel.proteinVertices{toBeChangedProtein(j)};
        eq = checkEqualityMatrices(vertices1,vertices2);
        % if j is being deleted anyway can't delete i
        
        if eq == 1 && ~contained(j,haveToBeDeletedProtein)
            if ~contained(i,haveToBeDeletedProtein)
                haveToBeDeletedProtein(index) = i;
                index = index + 1;
                break
                % break out of this loop now we know i is being deleted
                % anyway
            end
        end
    end
end

% same for LPS

toBeChangedLPS = setdiff(toBeChangedLPS,haveToBeDeletedLPS);

index = numel(haveToBeDeletedLPS)+1;
for i=1:numel(toBeChangedLPS)
    vertices1 = newModel.lpsVertices{toBeChangedLPS(i)};
    for j=i+1:numel(toBeChangedLPS)
        vertices2 = newModel.lpsVertices{toBeChangedLPS(j)};
        eq = checkEqualityMatrices(vertices1,vertices2);
        % if j is being deleted anyway can't delete i
        if eq == 1 && ~contained(j,haveToBeDeletedLPS)
            if ~contained(i,haveToBeDeletedLPS)
                haveToBeDeletedLPS(index) = i;
                index = index + 1;
                break
                % break out of this loop now we know i is being deleted
                % anyway
            end
        end
    end
end

% now we can delete everything in haveToBeDeleted and we should have no
% copies

haveToBeDeletedProtein;

newModel.proteinVertices;

newModel.proteinVertices(haveToBeDeletedProtein) = [];
newModel.lpsVertices(haveToBeDeletedLPS) = [];

newModel.proteinVertices;

% finally we add in the new regions

index = numel(newModel.proteinVertices)+1;

for i=1:numel(newVerticesProtein)
    newModel.proteinVertices{index} = newVerticesProtein{i};
    index = index + 1;
end

index = numel(newModel.lpsVertices)+1;
for i=1:numel(newVerticesLPS)
    newModel.lpsVertices{index} = newVerticesLPS{i};
    index = index + 1;
end

disp('BEFORE')
model.proteinVertices
disp('RESOLVED')
newModel.proteinVertices

end

% this tells us if a value is contained within an array and doesn't have
% issues with empty arrays
function result = contained(element,array)
    
    if numel(array) == 0
        result = 0;
        return
    end
    result = ismember(array,element);
end