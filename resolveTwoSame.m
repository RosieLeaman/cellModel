% this takes in two sets of vertices marking boundaries of two polygons of 
% the same type that are too close

% This code ASSUMES that the two polygons are only close in ONE PLACE. It
% WILL break if they are close in two places. It throws an error if this
% occurs

function newVertices = resolveTwoSame(vertices1,vertices2,plotYes)

figure; hold on;
plot(vertices1(:,1),vertices1(:,2),'x-')
plot(vertices2(:,1),vertices2(:,2),'x-')

vertices1Copy = vertices1;
vertices2Copy = vertices2;

% all that we know about these polygons is that they are close in at least
% one place. We need to find all the vertices where they are close.
% first we have to find the problem vertices. We store their INDICES.

problemVertices1 = [];
problemVertices2 = [];

index1 = 1;
index2 = 1;

for i=1:size(vertices1,1)
    x = vertices1(i,:);
    for j=1:size(vertices2,1)
        y = vertices2(j,:);
        if findDist(x,y) < 1
            disp('found one')
            % if they are too close, add the vertices as a problem
            % need to check whether this vertex has already been added and
            % not add it again
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

problemVertices1
problemVertices2

% sort the two arrays

problemVertices1 = sort(problemVertices1);
problemVertices2 = sort(problemVertices2);

% there are two different options.  Either the polygons intersect, or they
% are approaching each other. These can be differentiated by whether there
% is a large jump in the problem indices

% check for errors
test = diff(problemVertices1)
numBigJumps = nnz(find(test>1))
if numBigJumps == 0 || (numBigJumps==1 && problemVertices1(1) == 1)
    disp('approaching')
    
    polygon1Inside = inpolygon(vertices1(:,1),vertices1(:,2),vertices2(:,1),vertices2(:,2));
    
    sum(polygon1Inside)
    
    if sum(polygon1Inside) == 0
        % no intersections, approaching not inside.
        
        % make the polygons larger by 1. Fit a spline, obtain the normal,
        % and move all the points in vertices1Copy and vertices2Copy out by
        % 1*normal. Then the new polygons will intersect, and the same code
        % as below can be used with a few modifications
        
        [splineX1,splineY1] = fitSpline(vertices1);
        [~,normals1] = findTangentFromSplines(splineX1.breaks,splineX1,splineY1,1);
        [splineX2,splineY2] = fitSpline(vertices2);
        [~,normals2] = findTangentFromSplines(splineX2.breaks,splineX2,splineY2,1);
        for j=1:size(vertices1,1)
            vertices1Copy(j,:) = vertices1Copy(j,:)+normals1(j,:);
        end
        for j=1:size(vertices2,1)
            vertices2Copy(j,:) = vertices2Copy(j,:)+normals2(j,:);
        end
        
        % find the points in polygon 1 that are inside polygon 2. This function
        % returns a logical array
        polygon1Inside = inpolygon(vertices1Copy(:,1),vertices1Copy(:,2),vertices2Copy(:,1),vertices2Copy(:,2));
        polygon2Inside = inpolygon(vertices2Copy(:,1),vertices2Copy(:,2),vertices1Copy(:,1),vertices1Copy(:,2));

        figure; hold on;
        plot(vertices1Copy(:,1),vertices1Copy(:,2),'bx-')
        plot(vertices1Copy(polygon1Inside,1),vertices1Copy(polygon1Inside,2),'ko')
        plot(vertices2Copy(:,1),vertices2Copy(:,2),'rx-')
        plot(vertices2Copy(polygon2Inside,1),vertices2Copy(polygon2Inside,2),'ko')

        % find the points that are not an issue

        notProblemVertices1 = 1:size(vertices1Copy,1);
        notProblemVertices1 = notProblemVertices1(~polygon1Inside);
        notProblemVertices2 = 1:size(vertices2Copy,1);
        notProblemVertices2 = notProblemVertices2(~polygon2Inside);

        plot(vertices1Copy(notProblemVertices1,1),vertices1Copy(notProblemVertices1,2),'gs')
        plot(vertices2Copy(notProblemVertices2,1),vertices2Copy(notProblemVertices2,2),'gs')

        % connect the two
        % shuffle the indices so they are contiguous (in real space)
        test1 = diff(notProblemVertices1);
        test2 = diff(notProblemVertices2);

        N1 = find(test1 > 1);
        N2 = find(test2 > 1);

        if numel(N1) == 0 && numel(N2) > 0
            % polygon 1 is contiguous and 2 is not
            disp('case1')
            newVertices = [vertices2(notProblemVertices2(1:N2),:);vertices1(notProblemVertices1,:);vertices2(notProblemVertices2(N2+1:end),:)];
        elseif numel(N2) == 0 && numel(N1) > 0
            % opposite
            disp('case2')
            newVertices = [vertices1(notProblemVertices1(1:N1),:);vertices2(notProblemVertices2,:);vertices1(notProblemVertices1(N1+1:end),:)];
        else
            % neither is
            disp('case3')
            newVertices = [vertices1(notProblemVertices1(1:N1),:);vertices2(notProblemVertices2(N2+1:end),:);vertices2(notProblemVertices2(1:N2),:);vertices1(notProblemVertices1(N1+1:end),:)];
        end

        figure;
        plot(newVertices(:,1),newVertices(:,2),'bd-')
    end
    
elseif numBigJumps == 1
    disp('intersection')
    % find the points in polygon 1 that are inside polygon 2. This function
    % returns a logical array
    polygon1Inside = inpolygon(vertices1(:,1),vertices1(:,2),vertices2(:,1),vertices2(:,2));
    polygon2Inside = inpolygon(vertices2(:,1),vertices2(:,2),vertices1(:,1),vertices1(:,2));
    
    figure; hold on;
    plot(vertices1(:,1),vertices1(:,2),'bx-')
    plot(vertices1(polygon1Inside,1),vertices1(polygon1Inside,2),'ko')
    plot(vertices2(:,1),vertices2(:,2),'rx-')
    plot(vertices2(polygon2Inside,1),vertices2(polygon2Inside,2),'ko')
    
    % find the points that are not an issue
    
    notProblemVertices1 = 1:size(vertices1,1);
    notProblemVertices1 = notProblemVertices1(~polygon1Inside);
    notProblemVertices2 = 1:size(vertices2,1);
    notProblemVertices2 = notProblemVertices2(~polygon2Inside);
    
    plot(vertices1(notProblemVertices1,1),vertices1(notProblemVertices1,2),'gs')
    plot(vertices2(notProblemVertices2,1),vertices2(notProblemVertices2,2),'gs')
    
    % connect the two
    % shuffle the indices so they are contiguous (in real space)
    test1 = diff(notProblemVertices1);
    test2 = diff(notProblemVertices2);
    
    N1 = find(test1 > 1);
    N2 = find(test2 > 1);
    
    if numel(N1) == 0 && numel(N2) > 0
        % polygon 1 is contiguous and 2 is not
        disp('case1')
        newVertices = [vertices2(notProblemVertices2(1:N2),:);vertices1(notProblemVertices1,:);vertices2(notProblemVertices2(N2+1:end),:)];
    elseif numel(N2) == 0 && numel(N1) > 0
        % opposite
        disp('case2')
        newVertices = [vertices1(notProblemVertices1(1:N1),:);vertices2(notProblemVertices2,:);vertices1(notProblemVertices1(N1+1:end),:)];
    else
        % neither is
        disp('case3')
        newVertices = [vertices1(notProblemVertices1(1:N1),:);vertices2(notProblemVertices2(N2+1:end),:);vertices2(notProblemVertices2(1:N2),:);vertices1(notProblemVertices1(N1+1:end),:)];
    end
    
    figure;
    plot(newVertices(:,1),newVertices(:,2),'bd-')
end
% 
% % if we have no errors then let's continue
% % find not problematic vertices
% 
% notProblemVertices1 = setdiff(1:size(vertices1,1),problemVertices1);
% notProblemVertices2 = setdiff(1:size(vertices2,1),problemVertices2);
% 
% % store the sizes of both vertices
% 
% size1 = size(vertices1,1);
% size2 = size(vertices2,1);
% 
% % now we find the lowest index in problemVertices1
% 
% problemStartIndex1 = min(problemVertices1);
% problemEndIndex1 = max(problemVertices1);
% notProblemStartIndex1 = min(notProblemVertices1);
% notProblemEndIndex1 = max(notProblemVertices1);
% 
% if notProblemStartIndex1 == 1 && notProblemEndIndex1 == size1
%     startIndex1 = next(problemEndIndex1,size1);
%     endIndex1 = prev(problemStartIndex1,size1);
%     N1 = size1 - numel(problemVertices1);
%     notProblemIndices1 = nextN(startIndex1,size1,N1);
% else
%     notProblemIndices1 = notProblemStartIndex1:notProblemEndIndex1;
%     N1 = numel(notProblemVertices1);
% end
% 
% problemStartIndex2 = min(problemVertices2);
% problemEndIndex2 = max(problemVertices2);
% notProblemStartIndex2 = min(notProblemVertices2);
% notProblemEndIndex2 = max(notProblemVertices2);
% 
% if notProblemStartIndex2 == 1 && notProblemEndIndex2 == size1
%     startIndex2 = next(problemEndIndex2,size2);
%     endIndex2 = prev(problemStartIndex2,size2);
%     N2 = size2 - numel(problemVertices2);
%     notProblemIndices2 = nextN(startIndex2,size2,N2);
% else
%     notProblemIndices2 = notProblemStartIndex2:notProblemEndIndex2;
%     N2 = numel(notProblemVertices2);
% end
% 
% % we take the old vertices and whack them together in a new order
% 
% newVertices = [vertices1(notProblemIndices1,:);vertices2(notProblemIndices2,:)];
% 
% figure;
% plot(newVertices(:,1),newVertices(:,2),'x-')
% title('new')
% 
% assert(1==0)
% 
% % we shuffle them round slightly so that the holes are not at the end of
% % the polygon (makes the spline bad right where we want it good)
% 
% newIndices = shuffle(1:size(newVertices,1),floor(N1/2));
% newVertices = newVertices(newIndices,:);
% 
% % a figure that can be removed if we do not want it
% if plotYes == 1
%     figure;hold on;
%     for i=1:numel(problemVertices1)
%         plot(vertices1(problemVertices1(i),1),vertices1(problemVertices1(i),2),'b-o')
%     end
%     for i=1:numel(problemVertices2)
%         plot(vertices2(problemVertices2(i),1),vertices2(problemVertices2(i),2),'b-o')
%     end
%     plot(vertices1(:,1),vertices1(:,2),'r-x')
%     plot(vertices2(:,1),vertices2(:,2),'r-x')
%     plot(newVertices(:,1),newVertices(:,2),'k--x')
% end

% used to put a spline through it all, not sure if necessary.

% spline option 2 

% firstSplit = N1-floor(N1/2)+1;
% secondSplit = N1+N2-floor(N1/2)+1;
% sizeNew = size(newVertices,1);
% 
% xvals = newVertices(:,1);
% newXvals = [xvals(1:firstSplit-1);xvals(firstSplit)+0.5*(xvals(firstSplit+1)-xvals(firstSplit))];
% newXvals = [newXvals;xvals(firstSplit+2:secondSplit-1)];
% newXvals = [newXvals;xvals(secondSplit)+0.5*(xvals(secondSplit+1)-xvals(secondSplit))];
% newXvals = [newXvals;xvals(secondSplit+2:end)];
% 
% yvals = newVertices(:,2);
% newYvals = [yvals(1:firstSplit-1);yvals(firstSplit)+0.5*(yvals(firstSplit+1)-yvals(firstSplit))];
% newYvals = [newYvals;yvals(firstSplit+2:secondSplit-1)];
% newYvals = [newYvals;yvals(secondSplit)+0.5*(yvals(secondSplit+1)-yvals(secondSplit))];
% newYvals = [newYvals;yvals(secondSplit+2:end)];
% 
% % didn't divide by the number to make it between 0 to 1
% newt = 1:numel(newXvals);
% 
% % get the x and y splines separately
% ppvalX = spline(newt,newXvals);
% ppvalY = spline(newt,newYvals);
% 
% % we need some new t query points
% t2 = [1:(firstSplit-2),linspace(firstSplit-1,firstSplit+1,20)];
% t2 = [t2,firstSplit+2:secondSplit-3,linspace(secondSplit-2,secondSplit,20)];
% t2 = [t2,secondSplit+1:sizeNew];
% 
% ppx = ppval(ppvalX,t2);
% ppy = ppval(ppvalY,t2);
% 
% % new figure
% 
% if plotYes == 1
%     figure;
%     plot(newVertices(:,1),newVertices(:,2),'k--x')
%     hold on;
%     plot(ppx,ppy,'bd-')
% end
% 
% % ACTUALLY ASSIGN THE VERTICES
% 
% newVertices = [ppx',ppy'];

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