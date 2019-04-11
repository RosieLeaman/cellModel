% we have to split blobs when they become too stretched somewhere

function [newVertices1,newVertices2] = splitBlob(vertices,plotYes)

% first we must find the precise vertices that have issues

sizeVertices = size(vertices,1);

problemVertices = [];
index = 1;

for i=1:sizeVertices
    x = vertices(i,:);
    for j=i+1:sizeVertices
        y = vertices(j,:);
        if findDist(x,y) < 1
            if findDistRoundPolygon(i,j,vertices) > 2
                if numel(problemVertices) == 0
                    problemVertices(index) = i;
                    problemVertices(index+1) = j;
                    index = index + 2;
                else
                    if ~ismember(problemVertices,i)
                        problemVertices(index) = i;
                        index = index + 1;
                    end
                    if ~ismember(problemVertices,j)
                        problemVertices(index) = j;
                        index = index + 1;
                    end
                end
            end            
        end
    end
end

problemVertices = sort(problemVertices);

% find the not problemVertices
notProblemVertices = setdiff(1:size(vertices,1),problemVertices);

% find the contiguous bits

notProblemDiffs = diff(notProblemVertices);

indices = find(notProblemDiffs > 2);

if numel(indices) ~= 2
    error('The polygon was close to itself in three places somehow')
end

newIndices2 = notProblemVertices(indices(1)+1:indices(2));

newIndices1 = notProblemVertices(indices(2)+1:end);
newIndices1 = [newIndices1,notProblemVertices(1:indices(1))];

newVertices1 = vertices(newIndices1,:);
newVertices2 = vertices(newIndices2,:);

% we must shuffle the indices round so that splines aren't at the end where
% the hole is

newIndices1 = shuffle(1:numel(newIndices1),floor(numel(newIndices1)/2));
newIndices2 = shuffle(1:numel(newIndices2),floor(numel(newIndices2)/2));

newVertices1 = newVertices1(newIndices1,:);
newVertices2 = newVertices2(newIndices2,:);

if plotYes == 1
    figure; hold on;
    plot(vertices(:,1),vertices(:,2),'x-')

    plot(vertices(problemVertices,1),vertices(problemVertices,2),'ro')

    plot(newVertices1(:,1),newVertices1(:,2),'bs-')
    plot(newVertices2(:,1),newVertices2(:,2),'ks-')
end

% then we have to join up the holes

% join up vertices1

% first get new vertices with one in the middle of the hole

firstSplit = size(newVertices1,1)-floor(numel(newIndices1)/2)+1;

newXvals1 = newVertices1(1:(firstSplit-1),1);
newXvals1 = [newXvals1;newVertices1(firstSplit,1)+0.5*(newVertices1(firstSplit+1,1)-newVertices1(firstSplit,1))];
newXvals1 = [newXvals1;newVertices1(firstSplit+2:end,1)];
    
newYvals1 = newVertices1(1:(firstSplit-1),2);
newYvals1 = [newYvals1;newVertices1(firstSplit,2)+0.5*(newVertices1(firstSplit+1,2)-newVertices1(firstSplit,2))];
newYvals1 = [newYvals1;newVertices1(firstSplit+2:end,2)];

% didn't divide by the number to make it between 0 to 1
newt = 1:numel(newXvals1);

% get the x and y splines separately
ppvalX = spline(newt,newXvals1);
ppvalY = spline(newt,newYvals1);

% need new t query values

t2 = 1:firstSplit-2;
t2 = [t2,linspace(firstSplit-1,firstSplit+1,20)];
t2 = [t2,firstSplit+2:numel(newXvals1)];

ppx1 = ppval(ppvalX,t2);
ppy1 = ppval(ppvalY,t2);


% do the same for newVertices2

% first get new vertices with one in the middle of the hole

secondSplit = size(newVertices2,1)-floor(numel(newIndices2)/2)+1;

newXvals2 = newVertices2(1:(secondSplit-1),1);
newXvals2 = [newXvals2;newVertices2(secondSplit,1)+0.5*(newVertices2(secondSplit+1,1)-newVertices2(secondSplit,1))];
newXvals2 = [newXvals2;newVertices2(secondSplit+2:end,1)];
    
newYvals2 = newVertices2(1:(secondSplit-1),2);
newYvals2 = [newYvals2;newVertices2(secondSplit,2)+0.5*(newVertices2(secondSplit+1,2)-newVertices2(secondSplit,2))];
newYvals2 = [newYvals2;newVertices2(secondSplit+2:end,2)];

% didn't divide by the number to make it between 0 to 1
newt = 1:numel(newXvals2);

% get the x and y splines separately
ppvalX = spline(newt,newXvals2);
ppvalY = spline(newt,newYvals2);

% need new t query values

t2 = 1:secondSplit-2;
t2 = [t2,linspace(secondSplit-1,secondSplit+1,20)];
t2 = [t2,secondSplit+2:numel(newXvals2)];

ppx2 = ppval(ppvalX,t2);
ppy2 = ppval(ppvalY,t2);

if plotYes == 1
    figure; hold on;
    plot(vertices(:,1),vertices(:,2),'bx-')
    plot(ppx1,ppy1,'kx-')
    plot(ppx2,ppy2,'kx-')
end

% ACTUALLY ASSIGN THE NEW VERTICES
newVertices1 =  [ppx1',ppy1'];
newVertices2 = [ppx2',ppy2'];

end

% this is a helpful function
function previous = prev(x,maxi)
    if x > 1
        previous = x - 1;
    else
        previous = maxi;
    end  
end

function following = next(x,maxi)
    if x < maxi
        following = x + 1;
    else
        following = 1;
    end  
end

% produces the indices of the next N numbers going in a circle
function following = nextN(x,maxi,N)
    following = zeros(1,N);
    following(1) = x;
    for i=2:N
        following(i) = next(following(i-1),maxi);
    end

end

% shuffle the numbers in the 1xm list along by N, so that result(1) =
% list(N)
function result = shuffle(list,N)
    result = nextN(N,numel(list),numel(list));
end
