function mainModel

% model settings

polygonSides = 100;

membraneRadius = 1;

insRateProtein = 1;
insRateLPS = 1;

insRateBAM = 0.5;

time = 0;
dt = 0.01;
maxTime = 0.05;

materialAddedNewInsertion = insRateProtein*dt;

% set up the model

%BAMlocs = [[3,3];[3.3,3.3];[3.3,2.7]];
BAMlocs = [[3,3]];
LptDlocs = [];
proteinVertices = [];
lpsVertices = [];

% for all existing insertions add new protein material
for i=1:size(BAMlocs,1)
    vertices = findVerticesNewMaterial(BAMlocs(i,:),polygonSides,materialAddedNewInsertion);
    proteinVertices(:,:,i) = vertices;
end


% make a pretty plot
figure;hold on;
for poly = 1:size(proteinVertices,3)
    plot(proteinVertices(:,1,poly),proteinVertices(:,2,poly),'o')
end
xlim([0,4*pi*membraneRadius]);ylim([0,2*pi*membraneRadius])
plot(BAMlocs(:,1),BAMlocs(:,2),'xk','linewidth',2)

% while time is less than max time
while time < maxTime

    % see if we need to insert a new BAM
    
    newBamRandom = rand(1);
    
    % this should probably involve membrane area here somewhere
    if time < 0.01
    %if newBamRandom < insRateBAM*dt
    
        % if yes add a BAM to a new randomly chosen location
        
        newBAMLoc = [3,2.85];
%         newBAMLoc = rand(1,2); % gives uniform between 0,1
%         newBAMLoc(1) = newBAMLoc(1)*4*pi*membraneRadius;
%         newBAMLoc(2) = newBAMLoc(2)*2*pi*membraneRadius;
        
        newIndex = size(BAMlocs,1) + 1;
        
        BAMlocs(newIndex,:) = newBAMLoc;
        
        % add new material (protein)
        
        vertices = findVerticesNewMaterial(BAMlocs(newIndex,:),polygonSides,materialAddedNewInsertion);
        proteinVertices(:,:,newIndex) = vertices;
        
        plot(newBAMLoc(1),newBAMLoc(2),'xk','linewidth',2)
        plot(proteinVertices(:,1,newIndex),proteinVertices(:,2,newIndex),'ok')
        
    end
       
    % see if we need to add new lptD (or how many?)
    
        % if yes pick a BAM
        
        % add new material (lipid)
              
    % move polygon vertices
    % have to loop through all polygons, then all pairs of vertices
    for poly = 1:size(proteinVertices,3)
        for j=1:size(proteinVertices(:,:,1),1)

            flow = calcFlow(proteinVertices(j,:,poly),BAMlocs,insRateProtein,[],insRateLPS,membraneRadius);
            
            newPos = findNewVertexPosition(proteinVertices(j,:,poly),flow,dt);
            
            proteinVertices(j,:,poly) = newPos;
        end
    end
    
    % plot each time step
    %disp('plotting')
    for poly = 1:1
        plot(proteinVertices(:,1,poly),proteinVertices(:,2,poly),'o')
    end

    %fill(proteinVertices(:,1,1),proteinVertices(:,2,1),'r')
    
    
    % move insertion points (protein)
    
    for j = 1:size(BAMlocs,1)
        % find the flow from all points EXCEPT itself
        tempBAMlocs = BAMlocs;
        tempBAMlocs(j,:) = [];
        
        flow = calcFlow(BAMlocs(j,:),tempBAMlocs,insRateProtein,[],insRateLPS,membraneRadius);
            
        newPos = findNewVertexPosition(BAMlocs(j,:),flow,dt);
            
        BAMlocs(j,:) = newPos;
        
    end
    
    % move insertion points (lps)
    
    % resolve boundary issues
    
    
    % calculate anything you want to calculate
    
    
    % increase time
    
    time = time + dt;
    
end

for poly = 1:size(proteinVertices,3)
    plot(proteinVertices(:,1,poly),proteinVertices(:,2,poly),'x-')
end

%fill(proteinVertices(:,1,1),proteinVertices(:,2,1));
plot(BAMlocs(:,1),BAMlocs(:,2),'x','linewidth',2)
% legend('Initial','Step 1','Step 2','Step 3','Step 4','Step 5')
    