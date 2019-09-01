% this takes in two sets of vertices marking boundaries of two polygons of 
% the same type that are too close

% This code ASSUMES that the two polygons are only close in ONE PLACE. It
% WILL break if they are close in two places. It throws an error if this
% occurs

function newVertices = resolveTwoSame(vertices1,vertices2,plotYes)

if plotYes == 1
    figure; hold on;
    plot(vertices1(:,1),vertices1(:,2),'x-')
    plot(vertices2(:,1),vertices2(:,2),'x-')
end

vertices1Copy = vertices1;
vertices2Copy = vertices2;

% make the polygons larger by 1. Fit a spline, obtain the normal,
% and move all the points in vertices1Copy and vertices2Copy out by
% 1*normal. Then the new polygons will intersect, and the same code
% as below can be used with a few modifications

% see surfaceTensionPolygon for details about this normal finding
% code
numVertices1 = size(vertices1,1);
[splineX1,splineY1] = fitSpline([vertices1;vertices1;vertices1]);
[~,normals1] = findTangentFromSplines(splineX1.breaks,splineX1,splineY1,1);
normals1 = normals1(numVertices1+1:numVertices1*2,:);

numVertices2 = size(vertices2,1);
[splineX2,splineY2] = fitSpline([vertices2;vertices2;vertices2]);
[~,normals2] = findTangentFromSplines(splineX2.breaks,splineX2,splineY2,1);
normals2 = normals2(numVertices2+1:numVertices2*2,:);

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

if plotYes == 1
    figure; hold on;
    plot(vertices1Copy(:,1),vertices1Copy(:,2),'bx-')
    plot(vertices1Copy(polygon1Inside,1),vertices1Copy(polygon1Inside,2),'ko')
    plot(vertices2Copy(:,1),vertices2Copy(:,2),'rx-')
    plot(vertices2Copy(polygon2Inside,1),vertices2Copy(polygon2Inside,2),'ko')
end

% find the points that are not an issue

notProblemVertices1 = 1:size(vertices1Copy,1);
notProblemVertices1 = notProblemVertices1(~polygon1Inside);
notProblemVertices2 = 1:size(vertices2Copy,1);
notProblemVertices2 = notProblemVertices2(~polygon2Inside);

if plotYes == 1
    plot(vertices1Copy(notProblemVertices1,1),vertices1Copy(notProblemVertices1,2),'gs')
    plot(vertices2Copy(notProblemVertices2,1),vertices2Copy(notProblemVertices2,2),'gs')
end

% connect the two
% shuffle the indices so they are contiguous (in real space)
test1 = diff(notProblemVertices1);
test2 = diff(notProblemVertices2);

N1 = find(test1 > 1);
N2 = find(test2 > 1);

if numel(N1) == 0 && numel(N2) > 0
    % polygon 1 is contiguous and 2 is not
    %disp('case1')
    newVertices = [vertices2(notProblemVertices2(1:N2),:);vertices1(notProblemVertices1,:);vertices2(notProblemVertices2(N2+1:end),:)];
elseif numel(N2) == 0 && numel(N1) > 0
    % opposite
    %disp('case2')
    newVertices = [vertices1(notProblemVertices1(1:N1),:);vertices2(notProblemVertices2,:);vertices1(notProblemVertices1(N1+1:end),:)];
else
    % neither is
    %disp('case')
    newVertices = [vertices1(notProblemVertices1(1:N1),:);vertices2(notProblemVertices2(N2+1:end),:);vertices2(notProblemVertices2(1:N2),:);vertices1(notProblemVertices1(N1+1:end),:)];
end

if plotYes == 1
    figure;
    plot(newVertices(:,1),newVertices(:,2),'bd-')
end

% 
% % check for errors
% test = diff(problemVertices1)
% numBigJumps = nnz(find(test>1))
% if numBigJumps == 0 || (numBigJumps==1 && problemVertices1(1) == 1)
%     disp('approaching')
%     
%     polygon1Inside = inpolygon(vertices1(:,1),vertices1(:,2),vertices2(:,1),vertices2(:,2));
%     
%     sum(polygon1Inside)
%     
%     if sum(polygon1Inside) == 0
%         % no intersections, approaching not inside.
%         
%         % make the polygons larger by 1. Fit a spline, obtain the normal,
%         % and move all the points in vertices1Copy and vertices2Copy out by
%         % 1*normal. Then the new polygons will intersect, and the same code
%         % as below can be used with a few modifications
%         
%         % see surfaceTensionPolygon for details about this normal finding
%         % code
%         numVertices1 = size(vertices1,1);
%         [splineX1,splineY1] = fitSpline([vertices1;vertices1;vertices1]);
%         [~,normals1] = findTangentFromSplines(splineX1.breaks,splineX1,splineY1,1);
%         normals1 = normals1(numVertices1+1:numVertices1*2,:);
%         
%         numVertices2 = size(vertices2,1);
%         [splineX2,splineY2] = fitSpline([vertices2;vertices2;vertices2]);
%         [~,normals2] = findTangentFromSplines(splineX2.breaks,splineX2,splineY2,1);
%         normals2 = normals2(numVertices2+1:numVertices2*2,:);
%         
%         for j=1:size(vertices1,1)
%             vertices1Copy(j,:) = vertices1Copy(j,:)+normals1(j,:);
%         end
%         for j=1:size(vertices2,1)
%             vertices2Copy(j,:) = vertices2Copy(j,:)+normals2(j,:);
%         end
%         
%         % find the points in polygon 1 that are inside polygon 2. This function
%         % returns a logical array
%         polygon1Inside = inpolygon(vertices1Copy(:,1),vertices1Copy(:,2),vertices2Copy(:,1),vertices2Copy(:,2));
%         polygon2Inside = inpolygon(vertices2Copy(:,1),vertices2Copy(:,2),vertices1Copy(:,1),vertices1Copy(:,2));
% 
%         figure; hold on;
%         plot(vertices1Copy(:,1),vertices1Copy(:,2),'bx-')
%         plot(vertices1Copy(polygon1Inside,1),vertices1Copy(polygon1Inside,2),'ko')
%         plot(vertices2Copy(:,1),vertices2Copy(:,2),'rx-')
%         plot(vertices2Copy(polygon2Inside,1),vertices2Copy(polygon2Inside,2),'ko')
% 
%         % find the points that are not an issue
% 
%         notProblemVertices1 = 1:size(vertices1Copy,1);
%         notProblemVertices1 = notProblemVertices1(~polygon1Inside);
%         notProblemVertices2 = 1:size(vertices2Copy,1);
%         notProblemVertices2 = notProblemVertices2(~polygon2Inside);
% 
%         plot(vertices1Copy(notProblemVertices1,1),vertices1Copy(notProblemVertices1,2),'gs')
%         plot(vertices2Copy(notProblemVertices2,1),vertices2Copy(notProblemVertices2,2),'gs')
% 
%         % connect the two
%         % shuffle the indices so they are contiguous (in real space)
%         test1 = diff(notProblemVertices1);
%         test2 = diff(notProblemVertices2);
% 
%         N1 = find(test1 > 1);
%         N2 = find(test2 > 1);
% 
%         if numel(N1) == 0 && numel(N2) > 0
%             % polygon 1 is contiguous and 2 is not
%             disp('case1')
%             newVertices = [vertices2(notProblemVertices2(1:N2),:);vertices1(notProblemVertices1,:);vertices2(notProblemVertices2(N2+1:end),:)];
%         elseif numel(N2) == 0 && numel(N1) > 0
%             % opposite
%             disp('case2')
%             newVertices = [vertices1(notProblemVertices1(1:N1),:);vertices2(notProblemVertices2,:);vertices1(notProblemVertices1(N1+1:end),:)];
%         else
%             % neither is
%             disp('case3')
%             newVertices = [vertices1(notProblemVertices1(1:N1),:);vertices2(notProblemVertices2(N2+1:end),:);vertices2(notProblemVertices2(1:N2),:);vertices1(notProblemVertices1(N1+1:end),:)];
%         end
% 
%         figure;
%         plot(newVertices(:,1),newVertices(:,2),'bd-')
%     end
%     
% elseif numBigJumps == 1
%     disp('intersection')
%     % find the points in polygon 1 that are inside polygon 2. This function
%     % returns a logical array
%     polygon1Inside = inpolygon(vertices1(:,1),vertices1(:,2),vertices2(:,1),vertices2(:,2));
%     polygon2Inside = inpolygon(vertices2(:,1),vertices2(:,2),vertices1(:,1),vertices1(:,2));
%     
%     figure; hold on;
%     plot(vertices1(:,1),vertices1(:,2),'bx-')
%     plot(vertices1(polygon1Inside,1),vertices1(polygon1Inside,2),'ko')
%     plot(vertices2(:,1),vertices2(:,2),'rx-')
%     plot(vertices2(polygon2Inside,1),vertices2(polygon2Inside,2),'ko')
%     
%     % find the points that are not an issue
%     
%     notProblemVertices1 = 1:size(vertices1,1);
%     notProblemVertices1 = notProblemVertices1(~polygon1Inside);
%     notProblemVertices2 = 1:size(vertices2,1);
%     notProblemVertices2 = notProblemVertices2(~polygon2Inside);
%     
%     plot(vertices1(notProblemVertices1,1),vertices1(notProblemVertices1,2),'gs')
%     plot(vertices2(notProblemVertices2,1),vertices2(notProblemVertices2,2),'gs')
%     
%     % connect the two
%     % shuffle the indices so they are contiguous (in real space)
%     test1 = diff(notProblemVertices1);
%     test2 = diff(notProblemVertices2);
%     
%     N1 = find(test1 > 1);
%     N2 = find(test2 > 1);
%     
%     if numel(N1) == 0 && numel(N2) > 0
%         % polygon 1 is contiguous and 2 is not
%         disp('case1')
%         newVertices = [vertices2(notProblemVertices2(1:N2),:);vertices1(notProblemVertices1,:);vertices2(notProblemVertices2(N2+1:end),:)];
%     elseif numel(N2) == 0 && numel(N1) > 0
%         % opposite
%         disp('case2')
%         newVertices = [vertices1(notProblemVertices1(1:N1),:);vertices2(notProblemVertices2,:);vertices1(notProblemVertices1(N1+1:end),:)];
%     else
%         % neither is
%         disp('case3')
%         newVertices = [vertices1(notProblemVertices1(1:N1),:);vertices2(notProblemVertices2(N2+1:end),:);vertices2(notProblemVertices2(1:N2),:);vertices1(notProblemVertices1(N1+1:end),:)];
%     end
%     
%     figure;
%     plot(newVertices(:,1),newVertices(:,2),'bd-')
% end


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