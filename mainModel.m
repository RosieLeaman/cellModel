function mainModel(settings)

% model settings, set here if nargin < 1 (no inputs passed)

if nargin < 1
    settings.polygonSides = 200;

    settings.membraneCircumference = pi; % um
    settings.currentMaxLen = 2; %um
    settings.initialArea = settings.membraneCircumference*settings.currentMaxLen*2; % um^2

    settings.insRateProtein = 0.5; %um^2/s; estimate pi*0.0024*0.0024
    settings.insRateLPS = 0.5; % per LptD

    settings.insRateBAM = 0.5; %s^-1; estimate 7/60
    settings.insRateLptD = 0.5; % per BAM

    settings.growthRate = settings.insRateLPS*settings.insRateLptD*settings.insRateBAM + settings.insRateProtein*settings.insRateBAM; % units I think um^2/s^2
    settings.sqrtGrowthRate = sqrt(settings.growthRate); % units I think um/s

    settings.BAMsize = 0.01;

    settings.time = 0;
    settings.dt = 0.001; %s
    settings.maxTime = 2;

    settings.materialAddedNewInsertion = settings.insRateProtein*settings.dt;
end

polygonSides = settings.polygonSides;

membraneCircumference = settings.membraneCircumference; % um
currentMaxLen = settings.currentMaxLen; %um
initialArea = settings.initialArea; % um^2

insRateProtein = settings.insRateProtein; %um^2/s; estimate pi*0.0024*0.0024
insRateLPS = settings.insRateLPS; % per LptD

insRateBAM = settings.insRateBAM; %s^-1; estimate 7/60
insRateLptD = settings.insRateLptD; % per BAM

growthRate = settings.growthRate; % units I think um^2/s^2
sqrtGrowthRate = settings.sqrtGrowthRate; % units I think um/s

BAMsize = settings.BAMsize;

time = settings.time;
dt = settings.dt; %s
maxTime = settings.maxTime;

materialAddedNewInsertion = settings.materialAddedNewInsertion;

% set up the model (as a struct)

%BAMlocs = [[3,3];[3.3,3.3];[3.3,2.7]];
model.BAMlocs = [[0,0.5*pi]];
model.LptDlocs = [];
model.proteinVertices = [];
model.lpsVertices = [];

model.settings = settings;

% for all existing insertions add new protein material
for i=1:size(model.BAMlocs,1)
    vertices = findVerticesNewMaterial(model.BAMlocs(i,:),polygonSides,materialAddedNewInsertion);
    model.proteinVertices(:,:,i) = vertices;
end


% make a pretty plot
figure;hold on;

plot(model.BAMlocs(:,1),model.BAMlocs(:,2),'xk','linewidth',2)
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
        
        newBAMIndex = size(model.BAMlocs,1) + 1;
        
        model.BAMlocs(newBAMIndex,:) = newBAMloc;
        
        newBAMlocs = newBAMloc;
        
    end
       
    % see if we need to add new lptD 
    
    newLptDIndex = 1;
    
    for BAM = 1:size(model.BAMlocs,1)
        % for each BAM we roll a die and if it is less than some number we
        % add an LptD near that BAM
        
        a = rand(1);
        
        if a < insRateLptD*dt
            % we pick a random angle in [0,2pi)
            theta = rand(1)*2*pi;
            
            newLptDlocs(newLptDIndex,:) = [model.BAMlocs(BAM,1)+BAMsize*cos(theta),model.BAMlocs(BAM,2)+BAMsize*sin(theta)];
            
            newLptDIndex = newLptDIndex + 1;
        end
        
    end
    
    % add the new LptD to the LptD location list
    
    for lptD = 1:(newLptDIndex-1)
        newLptDLocsIndex = size(model.LptDlocs,1) + 1;

        model.LptDlocs(newLptDLocsIndex,:) = newLptDlocs(lptD,:);        
    end

    % move protein polygon vertices
    % have to loop through all polygons, then all pairs of vertices
    for poly = 1:size(model.proteinVertices,3)
        for j=1:size(model.proteinVertices(:,:,1),1)

            flow = calcFlow(model.proteinVertices(j,:,poly),model.BAMlocs,insRateProtein,model.LptDlocs,insRateLPS,membraneCircumference);
            
            newPos = findNewVertexPosition(model.proteinVertices(j,:,poly),flow,dt,membraneCircumference);
            
            model.proteinVertices(j,:,poly) = newPos;
        end
    end
    
    % move lps polygon vertices
    % have to loop through all polygons, then all pairs of vertices
    for poly = 1:size(model.lpsVertices,3)
        for j=1:size(model.lpsVertices(:,:,1),1)

            flow = calcFlow(model.lpsVertices(j,:,poly),model.BAMlocs,insRateProtein,model.LptDlocs,insRateLPS,membraneCircumference);
            
            newPos = findNewVertexPosition(model.lpsVertices(j,:,poly),flow,dt,membraneCircumference);
            
            model.lpsVertices(j,:,poly) = newPos;
        end
    end
    
    % move insertion points (protein)
    
    for j = 1:size(model.BAMlocs,1)
        % find the flow from all points EXCEPT itself
        tempBAMlocs = model.BAMlocs;
        tempBAMlocs(j,:) = [];
        
        flow = calcFlow(model.BAMlocs(j,:),tempBAMlocs,insRateProtein,model.LptDlocs,insRateLPS,membraneCircumference);
            
        newPos = findNewVertexPosition(model.BAMlocs(j,:),flow,dt,membraneCircumference);
            
        model.BAMlocs(j,:) = newPos;
        
    end
    
    % move insertion points (lps)
    
    for j = 1:size(model.LptDlocs,1)
        % find the flow from all points EXCEPT itself
        tempLptDlocs = model.LptDlocs;
        tempLptDlocs(j,:) = [];
        
        flow = calcFlow(model.LptDlocs(j,:),model.BAMlocs,insRateProtein,tempLptDlocs,insRateLPS,membraneCircumference);
            
        newPos = findNewVertexPosition(model.LptDlocs(j,:),flow,dt,membraneCircumference);
            
        model.LptDlocs(j,:) = newPos;
        
    end
    
    % add new material from insertions
    
    % add new material (protein)
    % check to see whether new BAMs were added or not
    if ~isempty(newBAMlocs)
        % loop through them (in case somehow we have more than one added)
        for j=1:size(newBAMlocs,1)
            vertices = findVerticesNewMaterial(newBAMlocs(j,:),polygonSides,materialAddedNewInsertion);
            model.proteinVertices(:,:,newBAMIndex) = vertices;

            plot(newBAMlocs(j,1),newBAMlocs(j,2),'xk','linewidth',2)
           
        end
    end
    
    % add new material (LPS)
    % check to see whether new LptDs were added or not
    if ~isempty(newLptDlocs)
        % loop through them (in case somehow we have more than one added)
        for j=1:size(newLptDlocs,1)
            vertices = findVerticesNewMaterial(newLptDlocs(j,:),polygonSides,materialAddedNewInsertion);
            
            newLptDIndex = size(model.lpsVertices,3) + 1;
            
            model.lpsVertices(:,:,newLptDIndex) = vertices;

            plot(newLptDlocs(j,1),newLptDlocs(j,2),'dk','linewidth',2)
           
        end
    end
        
    
    % resolve boundary issues
    
    
    % calculate anything you want to calculate
    
    
    % increase time
    
    time = time + dt;
    
end

for poly = 1:size(model.proteinVertices,3)
    plot(model.proteinVertices(:,1,poly),model.proteinVertices(:,2,poly),'x')
end

plot(model.BAMlocs(:,1),model.BAMlocs(:,2),'x','linewidth',2)


if numel(model.LptDlocs) > 0
    plot(model.LptDlocs(:,1),model.LptDlocs(:,2),'d','linewidth',2)
end

visualise(model);
save('results.mat','model')