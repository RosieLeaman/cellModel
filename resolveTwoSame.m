% this takes in two sets of vertices marking boundaries of two polygons of 
% the same type that are too close

% This code ASSUMES that the two polygons are only close in ONE PLACE. It
% WILL break if they are close in two places. It throws an error if this
% occurs

function newVertices = resolveTwoSame(vertices1,vertices2)

% first we have to find the problem vertices. We store their INDICES

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

% sort the two arrays

problemVertices1 = sort(problemVertices1);
problemVertices2 = sort(problemVertices2);

% check for errors
test = diff(problemVertices1);
if nnz(find(test>1)) > 1
    error('There were two polygons which were close at two points, not 1')
end

test = diff(problemVertices2);
if nnz(find(test>1)) > 1
    error('There were two polygons which were close at two points, not 1')
end

% if we have no errors then let's continue
% find not problematic vertices

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

newVertices = [vertices1(notProblemIndices1,:);vertices2(notProblemIndices2,:)];

% we shuffle them round slightly so that the holes are not at the end of
% the polygon (makes the spline bad right where we want it good)

newIndices = shuffle(1:size(newVertices,1),floor(N1/2));
newVertices = newVertices(newIndices,:);

% a figure that can be removed if we do not want it
figure;hold on;
plot(vertices1(:,1),vertices1(:,2),'r-x')
plot(vertices2(:,1),vertices2(:,2),'r-x')

for i=1:numel(problemVertices1)
    plot(vertices1(problemVertices1(i),1),vertices1(problemVertices1(i),2),'b-o')
end
for i=1:numel(problemVertices2)
    plot(vertices2(problemVertices2(i),1),vertices2(problemVertices2(i),2),'b-o')
end

plot(newVertices(:,1),newVertices(:,2),'k--x')

% now we want to whack a spline through all that.

% new figure

figure;
plot(newVertices(:,1),newVertices(:,2),'k--x')

% first get the distances apart from each other
sizeNew = size(newVertices,1);
t = 1:sizeNew;

% normalise to [0,1]
t = t./(max(t));

% get the x and y splines separately
ppvalX = spline(t,newVertices(:,1));
ppvalY = spline(t,newVertices(:,2));

% we need some new t query points
firstSplit = N1-floor(N1/2)+1;
secondSplit = N1+N2-floor(N1/2)+1;
t2 = [1:(firstSplit-1),linspace(firstSplit,firstSplit+1,20)];
t2 = [t2,firstSplit+2:secondSplit-1,linspace(secondSplit,secondSplit+1,20)];
t2 = [t2,secondSplit+2:sizeNew];

t2 = t2./(max(t2));
 

ppx = ppval(ppvalX,t2);
ppy = ppval(ppvalY,t2);

hold on;
plot(ppx,ppy,'ro-')

% spline option 2

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

hold on;

% we need some new t query points
t2 = [1:(firstSplit-2),linspace(firstSplit-1,firstSplit+1,20)];
t2 = [t2,firstSplit+2:secondSplit-3,linspace(secondSplit-2,secondSplit,20)];
t2 = [t2,secondSplit+1:sizeNew];

 

ppx = ppval(ppvalX,t2);
ppy = ppval(ppvalY,t2);

hold on;
plot(ppx,ppy,'bd-')

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