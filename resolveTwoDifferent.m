% this resolves the issues with a region interior to another reaching the
% edge

function [newVertices,insideIndex] = resolveTwoDifferent(vertices1,vertices2,plotYes)

% first we have to work out which one is inside the other
% this could not be done using pointInPolygon. It can be done using mins
% and maxes

max1 = max(vertices1(:,1));
min1 = min(vertices1(:,1));
max2 = max(vertices2(:,1));
min2 = min(vertices2(:,1));

% test if vertices2 are inside vertices 1

if min2 > min1 && max2 < max1
    inside = 1;
else
    inside = 0;
end

% we will arrange it so that vertices2 is inside and vertices1 is outside

if inside ~= 1
    % if this is not currently the case we have to swap them
    temp = vertices1;
    vertices1 = vertices2;
    vertices2 = temp;
    insideIndex = 1;
    disp('deleting protein')
else
    insideIndex = 2;
    disp('deleting lps')
end

% vertices1 will be edited to add the points from vertices2 and vertices2
% will be deleted

% the above should actually be removed and put elsewhere and we should make
% the ASSUMPTION that the inputs are given to us ordered with vertices1
% inside and vertices2 outside

% why should this be removed outside?

% now we can merrily go along on our usual approach. Same as resolveTwoSame
% we first find the issues

problemVertices1 = [];
problemVertices2 = [];

index1 = 1;
index2 = 1;

for i=1:size(vertices1,1)
    x = vertices1(i,:);
    for j=1:size(vertices2,1)
        y = vertices2(j,:);
        if findDist(x,y) < 1
            if numel(problemVertices1) == 0
                problemVertices1(index1) = i;
                index1 = index1 + 1;
            elseif ~ismember(problemVertices1,i)
                problemVertices1(index1) = i;
                index1 = index1 + 1;
            end
            
            if numel(problemVertices2) == 0
                problemVertices2(index2) = j;
                index2 = index2 + 1;
            elseif ~ismember(problemVertices2,j)
                problemVertices2(index2) = j;
                index2 = index2 + 1;
            end
        end
    end
end

% have an escape here in case something has changed and there aren't any
% issues any more

if numel(problemVertices1) == 0
    error('Something has gone wrong, these are not actually close')
end

% sort the two arrays

problemVertices1 = sort(problemVertices1);
problemVertices2 = sort(problemVertices2);

% check that there is only one place where they are close

test = diff(problemVertices1);
if nnz(find(test>1)) > 1
    error('These polygons were close at two points, not 1')
end

test = diff(problemVertices2);
if nnz(find(test>1)) > 1
    error('These polygons were close at two points, not 1')
end

notProblemVertices1 = setdiff(1:size(vertices1,1),problemVertices1);
notProblemVertices2 = setdiff(1:size(vertices2,1),problemVertices2);

% store the sizes of both vertices

size1 = size(vertices1,1);
size2 = size(vertices2,1);

% now we find the lowest index in problemVertices1

problemStartIndex1 = min(problemVertices1);
problemEndIndex1 = max(problemVertices1);
notProblemStartIndex1 = min(notProblemVertices1);
notProblemEndIndex1 = max(notProblemVertices1);

if notProblemStartIndex1 == 1 && notProblemEndIndex1 == size1
    startIndex1 = next(problemEndIndex1,size1);
    endIndex1 = prev(problemStartIndex1,size1);
    N1 = size1 - numel(problemVertices1);
    notProblemIndices1 = nextN(startIndex1,size1,N1);
else
    notProblemIndices1 = notProblemStartIndex1:notProblemEndIndex1;
    N1 = numel(notProblemVertices1);
end

problemStartIndex2 = min(problemVertices2);
problemEndIndex2 = max(problemVertices2);
notProblemStartIndex2 = min(notProblemVertices2);
notProblemEndIndex2 = max(notProblemVertices2);

if notProblemStartIndex2 == 1 && notProblemEndIndex2 == size1
    startIndex2 = next(problemEndIndex2,size2);
    endIndex2 = prev(problemStartIndex2,size2);
    N2 = size2 - numel(problemVertices2);
    notProblemIndices2 = nextN(startIndex2,size2,N2);
else
    notProblemIndices2 = notProblemStartIndex2:notProblemEndIndex2;
    N2 = numel(notProblemVertices2);
end

% we take the old vertices and whack them together in a new order

% in this case we have to reverse the order of the indices for vertices2
% note this is the ONLY difference between this and resolveTwoSame
notProblemIndices2Rev = fliplr(notProblemIndices2);

newVertices = [vertices1(notProblemIndices1,:);vertices2(notProblemIndices2Rev,:)];

% we shuffle them round slightly so that the holes are not at the end of
% the polygon (makes the spline bad right where we want it good)

newIndices = shuffle(1:size(newVertices,1),floor(N1/2));
newVertices = newVertices(newIndices,:);

for i=1:numel(problemVertices1)
    plot(vertices1(problemVertices1(i),1),vertices1(problemVertices1(i),2),'b-o')
end
for i=1:numel(problemVertices2)
    plot(vertices2(problemVertices2(i),1),vertices2(problemVertices2(i),2),'b-o')
end

% a figure that can be removed if we do not want it
if plotYes == 1
    figure;hold on;
    plot(vertices1(:,1),vertices1(:,2),'r-x')
    plot(vertices2(:,1),vertices2(:,2),'r-x')
    plot(newVertices(:,1),newVertices(:,2),'k--x')
end


% now we want to whack a spline through all that.


% % first get the distances apart from each other
% sizeNew = size(newVertices,1);
% t = 1:sizeNew;
% 
% % normalise to [0,1]
% t = t./(max(t));
% 
% % get the x and y splines separately
% ppvalX = spline(t,newVertices(:,1));
% ppvalY = spline(t,newVertices(:,2));
% 
% % we need some new t query points
% firstSplit = N1-floor(N1/2)+1;
% secondSplit = N1+N2-floor(N1/2)+1;
% t2 = [1:(firstSplit-1),linspace(firstSplit,firstSplit+1,20)];
% t2 = [t2,firstSplit+2:secondSplit-1,linspace(secondSplit,secondSplit+1,20)];
% t2 = [t2,secondSplit+2:sizeNew];
% 
% t2 = t2./(max(t2));
%  
% 
% ppx = ppval(ppvalX,t2);
% ppy = ppval(ppvalY,t2);
% 
% hold on;
% plot(ppx,ppy,'ro-')

% spline option 2 THIS OPTION IS USED

firstSplit = N1-floor(N1/2)+1;
secondSplit = N1+N2-floor(N1/2)+1;
sizeNew = size(newVertices,1);

xvals = newVertices(:,1);
newXvals = [xvals(1:firstSplit-1);xvals(firstSplit)+0.5*(xvals(firstSplit+1)-xvals(firstSplit))];
newXvals = [newXvals;xvals(firstSplit+2:secondSplit-1)];
newXvals = [newXvals;xvals(secondSplit)+0.5*(xvals(secondSplit+1)-xvals(secondSplit))];
newXvals = [newXvals;xvals(secondSplit+2:end)];

yvals = newVertices(:,2);
newYvals = [yvals(1:firstSplit-1);yvals(firstSplit)+0.5*(yvals(firstSplit+1)-yvals(firstSplit))];
newYvals = [newYvals;yvals(firstSplit+2:secondSplit-1)];
newYvals = [newYvals;yvals(secondSplit)+0.5*(yvals(secondSplit+1)-yvals(secondSplit))];
newYvals = [newYvals;yvals(secondSplit+2:end)];

% didn't divide by the number to make it between 0 to 1
newt = 1:numel(newXvals);

% get the x and y splines separately
ppvalX = spline(newt,newXvals);
ppvalY = spline(newt,newYvals);

% we need some new t query points
t2 = [1:(firstSplit-2),linspace(firstSplit-1,firstSplit+1,20)];
t2 = [t2,firstSplit+2:secondSplit-3,linspace(secondSplit-2,secondSplit,20)];
t2 = [t2,secondSplit+1:sizeNew];

 
ppx = ppval(ppvalX,t2);
ppy = ppval(ppvalY,t2);

% new figure
if plotYes == 1
    figure;
    plot(newVertices(:,1),newVertices(:,2),'k--x')
    hold on;
    plot(ppx,ppy,'bd-')
end

% ACTUALLY ASSIGN THE VERTICES

newVertices = [ppx',ppy'];

size(newVertices)

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