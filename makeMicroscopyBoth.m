function [newGridLPS,newGridProtein] = makeMicroscopyBoth(lpsVertices,proteinVertices,xMax,yMax)

% plot the thing prettily so we know what it looks like
visualisePretty(proteinVertices,lpsVertices)

% order the LPS vertices in groups which have vertices in rolling chunks

% make a grid

gridLPS = zeros(yMax+401,2*xMax+1);
gridProtein = zeros(yMax+401,2*xMax+1);
xVals = -xMax:xMax;
yVals = -400:yMax;

limits = linspace(-xMax,xMax,10);
currentLimitIndex = 2;
currentLimit = limits(currentLimitIndex);

% make a mini lps vertices
miniLpsVertices = {};
index = 1;
for k=1:numel(lpsVertices)
    if sum(lpsVertices{k}(:,1) < limits(currentLimitIndex)) && sum(lpsVertices{k}(:,1) > limits(currentLimitIndex-1))
        miniLpsVertices{index} = lpsVertices{k};
        index = index + 1;
    end
end
miniProteinVertices = {};
index = 1;
for k=1:numel(proteinVertices)
    if sum(proteinVertices{k}(:,1) < limits(currentLimitIndex)) && sum(proteinVertices{k}(:,1) > limits(currentLimitIndex-1))
        miniProteinVertices{index} = proteinVertices{k};
        index = index + 1;
    end
end

gridSizeX = size(gridLPS,2);
gridSizeY = size(gridLPS,1);

for i=1:gridSizeX
    % go along x
    
    % check whether we need to change our current limit
    if xVals(i) > currentLimit
        % expand limit
        currentLimitIndex = currentLimitIndex + 1;
        currentLimit = limits(currentLimitIndex);
        
        % recalculate nearby lps polygons
        % make a mini lps vertices
        miniLpsVertices = {};
        index = 1;
        for k=1:numel(lpsVertices)
            if sum(lpsVertices{k}(:,1) < limits(currentLimitIndex)) && sum(lpsVertices{k}(:,1) > limits(currentLimitIndex-1))
                miniLpsVertices{index} = lpsVertices{k};
                index = index + 1;
            end
        end
        miniProteinVertices = {};
        index = 1;
        for k=1:numel(proteinVertices)
            if sum(proteinVertices{k}(:,1) < limits(currentLimitIndex)) && sum(proteinVertices{k}(:,1) > limits(currentLimitIndex-1))
                miniProteinVertices{index} = proteinVertices{k};
                index = index + 1;
            end
        end
    end
    
    % go through each polygon find if in LPS
    columnLPS = zeros(gridSizeY,1);
    for k=1:numel(miniLpsVertices)
        insideLPS = inpolygon(xVals(i)*ones(numel(yVals),1),yVals,miniLpsVertices{k}(:,1),miniLpsVertices{k}(:,2));
    
        columnLPS = columnLPS + double(insideLPS);
    end
    
    % find if in protein
    columnProtein = zeros(gridSizeY,1);
    for k=1:numel(miniProteinVertices)
        insideProtein = inpolygon(xVals(i)*ones(numel(yVals),1),yVals,miniProteinVertices{k}(:,1),miniProteinVertices{k}(:,2));

        columnProtein = columnProtein + double(insideProtein);
    end
    
    actualColumnLPS = columnLPS;
    
    actualColumnLPS(columnLPS < columnProtein) = 0;
    actualColumnLPS(columnLPS >= columnProtein & columnLPS > 0) = 1;
    
    gridLPS(:,i) = actualColumnLPS;
    
    actualColumnProtein = columnProtein;
    
    actualColumnProtein(columnProtein <= columnLPS) = 0;
    actualColumnProtein(columnProtein > columnLPS) = 1;

    gridProtein(:,i) = actualColumnProtein;
        
end

% figure;
% spy(grid)

pxSize = 97;

gridLPS = imgaussfilt(gridLPS,pxSize);
gridProtein = imgaussfilt(gridProtein,pxSize);

% figure;
% imagesc(gridProtein)

% figure;
% imagesc(grid);

% bin the grid smaller

newGridSizeX = floor(gridSizeX/pxSize);
newGridSizeY = floor(gridSizeY/pxSize);

newGridLPS = zeros(newGridSizeY,newGridSizeX);

for i=1:newGridSizeX
    for j=1:newGridSizeY
        smolGrid = gridLPS(1+pxSize*(j-1):pxSize*j,1+pxSize*(i-1):pxSize*i);
        newGridLPS(j,i) = mean(smolGrid(:));   
    end
end

% flip the y axis so that it matches the original figure
for i=1:newGridSizeX
    column = newGridLPS(:,i);
    column = column(newGridSizeY:-1:1);
    newGridLPS(:,i) = column;
end


newGridProtein = zeros(newGridSizeY,newGridSizeX);

for i=1:newGridSizeX
    for j=1:newGridSizeY
        smolGrid = gridProtein(1+pxSize*(j-1):pxSize*j,1+pxSize*(i-1):pxSize*i);
        newGridProtein(j,i) = mean(smolGrid(:));   
    end
end

% flip the y axis so that it matches the original figure
for i=1:newGridSizeX
    column = newGridProtein(:,i);
    column = column(newGridSizeY:-1:1);
    newGridProtein(:,i) = column;
end

% figure;
% imagesc(newGridProtein)

figure;
imshow(newGridLPS,[min(newGridLPS(:)),max(newGridLPS(:))]);
xticklabels({''});
yticklabels({''});

figure;
imshow(newGridProtein,[min(newGridProtein(:)),max(newGridProtein(:))]);
xticklabels({''});
yticklabels({''});

figure;
plot(mean(newGridLPS),'linewidth',2);
axis tight
xlabel('Distance along long axis')
ylabel('Average fluorescence intensity')
set(gca,'FontSize',24)
grid on
xticklabels({''})

figure;
plot(mean(newGridProtein),'linewidth',2);
axis tight
xlabel('Distance along long axis')
ylabel('Average fluorescence intensity')
set(gca,'FontSize',24)
grid on
xticklabels({''})

end
