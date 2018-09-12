function mainModel

% model settings

polygonSides = 200;

membraneCircumference = pi; % um
currentMaxLen = 2; %um
initialArea = membraneCircumference*currentMaxLen*2; % um^2

insRateProtein = 0.5; %um^2/s; estimate pi*0.0024*0.0024
insRateLPS = 0.5; % per LptD

insRateBAM = 0.1; %s^-1; estimate 7/60
insRateLptD = 0.005; % per BAM

growthRate = insRateLPS*insRateLptD*insRateBAM + insRateProtein*insRateBAM; % units I think um^2/s^2
sqrtGrowthRate = sqrt(growthRate); % units I think um/s

BAMsize = 0.01;

time = 0;
dt = 0.01; %s
maxTime = 4;

materialAddedNewInsertion = insRateProtein*dt;

% set up the model

%BAMlocs = [[3,3];[3.3,3.3];[3.3,2.7]];
BAMlocs = [[0,0.5*pi]];
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

plot(BAMlocs(:,1),BAMlocs(:,2),'xk','linewidth',2)
plot([currentMaxLen,currentMaxLen],[0,membraneCircumference],'k-')
plot([-currentMaxLen,-currentMaxLen],[0,membraneCircumference],'k-')

% while time is less than max time
while time < maxTime

    % we need to maintain a list of newly inserted insertion points, or at
    % least their indices
    % there can only be at most one new bam location, and this is already
    % noted
    
    newBAMlocs = [];
    newLptDlocs = [];
    
    % see if we need to insert a new BAM
    
    newBamRandom = rand(1);
    
    % work out current area
    
    currentArea = initialArea*exp(sqrtGrowthRate*time);
    
    if newBamRandom < insRateBAM*dt*currentArea
        disp(time)
        % if yes add a BAM to a new randomly chosen location
        
        newBAMloc = rand(1,2); % gives two uniform random numbers
        
        % adjust the x one to be uniform between -currentMaxLen and
        % currentMaxLen
        % and y one to be uniform between 0 and membraneCircumference
        newBAMloc(1) = newBAMloc(1)*2*currentMaxLen-currentMaxLen;
        newBAMloc(2) = newBAMloc(2)*membraneCircumference;
        
        newBAMIndex = size(BAMlocs,1) + 1;
        
        BAMlocs(newBAMIndex,:) = newBAMloc;
        
        newBAMlocs = newBAMloc;
        
    end
       
    % see if we need to add new lptD 
    
    newLptDIndex = 1;
    
    for BAM = 1:size(BAMlocs,1)
        % for each BAM we roll a die and if it is less than some number we
        % add an LptD near that BAM
        
        a = rand(1);
        
        if a < insRateLptD
            % we pick a random angle in [0,2pi)
            theta = rand(1)*2*pi;
            
            newLptDlocs(newLptDIndex,:) = [BAMlocs(BAM,1)+BAMsize*cos(theta),BAMlocs(BAM,2)+BAMsize*sin(theta)];
            
            newLptDIndex = newLptDIndex + 1;
        end
        
    end
    
    % add the new LptD to the LptD location list
    
    for lptD = 1:(newLptDIndex-1)
        newLptDLocsIndex = size(LptDlocs,1) + 1;

        LptDlocs(newLptDLocsIndex,:) = newLptDlocs(lptD,:);        
    end

    % move protein polygon vertices
    % have to loop through all polygons, then all pairs of vertices
    for poly = 1:size(proteinVertices,3)
        for j=1:size(proteinVertices(:,:,1),1)

            flow = calcFlow(proteinVertices(j,:,poly),BAMlocs,insRateProtein,LptDlocs,insRateLPS,membraneCircumference);
            
            newPos = findNewVertexPosition(proteinVertices(j,:,poly),flow,dt,membraneCircumference);
            
            proteinVertices(j,:,poly) = newPos;
        end
    end
    
    % move lps polygon vertices
    % have to loop through all polygons, then all pairs of vertices
    for poly = 1:size(lpsVertices,3)
        for j=1:size(lpsVertices(:,:,1),1)

            flow = calcFlow(lpsVertices(j,:,poly),BAMlocs,insRateProtein,LptDlocs,insRateLPS,membraneCircumference);
            
            newPos = findNewVertexPosition(lpsVertices(j,:,poly),flow,dt,membraneCircumference);
            
            lpsVertices(j,:,poly) = newPos;
        end
    end
    
    % move insertion points (protein)
    
    for j = 1:size(BAMlocs,1)
        % find the flow from all points EXCEPT itself
        tempBAMlocs = BAMlocs;
        tempBAMlocs(j,:) = [];
        
        flow = calcFlow(BAMlocs(j,:),tempBAMlocs,insRateProtein,LptDlocs,insRateLPS,membraneCircumference);
            
        newPos = findNewVertexPosition(BAMlocs(j,:),flow,dt,membraneCircumference);
            
        BAMlocs(j,:) = newPos;
        
    end
    
    % move insertion points (lps)
    
    for j = 1:size(LptDlocs,1)
        % find the flow from all points EXCEPT itself
        tempLptDlocs = LptDlocs;
        tempLptDlocs(j,:) = [];
        
        flow = calcFlow(LptDlocs(j,:),BAMlocs,insRateProtein,tempLptDlocs,insRateLPS,membraneCircumference);
            
        newPos = findNewVertexPosition(LptDlocs(j,:),flow,dt,membraneCircumference);
            
        LptDlocs(j,:) = newPos;
        
    end
    
    % add new material from insertions
    
    % add new material (protein)
    % check to see whether new BAMs were added or not
    if ~isempty(newBAMlocs)
        % loop through them (in case somehow we have more than one added)
        for j=1:size(newBAMlocs,1)
            vertices = findVerticesNewMaterial(newBAMlocs(j,:),polygonSides,materialAddedNewInsertion);
            proteinVertices(:,:,newBAMIndex) = vertices;

            plot(newBAMlocs(j,1),newBAMlocs(j,2),'xk','linewidth',2)
           
        end
    end
    
    % add new material (LPS)
    % check to see whether new LptDs were added or not
    if ~isempty(newLptDlocs)
        % loop through them (in case somehow we have more than one added)
        for j=1:size(newLptDlocs,1)
            vertices = findVerticesNewMaterial(newLptDlocs(j,:),polygonSides,materialAddedNewInsertion);
            
            newLptDIndex = size(lpsVertices,3) + 1;
            
            lpsVertices(:,:,newLptDIndex) = vertices;

            plot(newLptDlocs(j,1),newLptDlocs(j,2),'dk','linewidth',2)
           
        end
    end
        
    
    % resolve boundary issues
    
    
    % calculate anything you want to calculate
    
    
    % increase time
    
    time = time + dt;
    
end

for poly = 1:size(proteinVertices,3)
    plot(proteinVertices(:,1,poly),proteinVertices(:,2,poly),'x')
end

plot(BAMlocs(:,1),BAMlocs(:,2),'x','linewidth',2)


if numel(LptDlocs) > 0
    plot(LptDlocs(:,1),LptDlocs(:,2),'d','linewidth',2)
end

figure; hold on;
for poly = 1:size(proteinVertices,3)
    h = fill(proteinVertices(:,1,poly),proteinVertices(:,2,poly),'b');
    set(h,'facealpha',.5)
end
if numel(lpsVertices) > 0
    for poly = 1:size(lpsVertices,3)
        h = fill(lpsVertices(:,1,poly),lpsVertices(:,2,poly),'r');
        set(h,'facealpha',.5)
    end
    plot(LptDlocs(:,1),LptDlocs(:,2),'d','linewidth',2)
end
plot(BAMlocs(:,1),BAMlocs(:,2),'x','linewidth',2)

   
