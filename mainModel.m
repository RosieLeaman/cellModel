% This is the main code for the model. It takes given positions and runs
% the code for so many iterations

% INPUTS:
% plotYes; 0 or 1 depending on whether we want to plot figures or not
% settings; a struct containing data about the model variables. Required
% fields are
%   polygonSides; number of sides to use to approximate polygons
%   membraneCircumference; the circumference of the cylindrical part of the
%       membrane
%   currentMaxLen; the initial length of the membrane (determines positions
%       of new BAM)
%   insRateProtein; insertion rate of protein from BAM
%   insRateLPS; insertion rate of LPS from LptD
%   insRateBAM; insertion rate of BAM complexes
%   insRateLptD; insertion rate of LptD from BAM
%   BAMsize; size of a BAM complex (radius)
%   time; initial time (should basically always be zero)
%   dt; time step size 
%   maxTime; time at which simulation will end (number of steps is
%       (maxTime-time)/dt)
%   surfaceTensionFlag; 0 or 1 depending if you want surface tension to be
%       included in the model
%   saveLocation; folder location to save any information at
%   for more details see the code below
%
% initPositions; a struct that determines the initial positions of BAM/LptD
% required fields are:
%   BAMlocs; nx2 matrix of [x,y] row vectors containing locations of BAM complexes
%   LptDlocs; mx2 matrix of [x,y] row vectors containing locations of LptD complexes
%   proteinVertices; nxkx2 matrix where n is number of BAM complexes, k is
%       number of vertices for each polygon and each entry is a row vector
%       [x,y] with position of each vertex
%   lpsVertices; mxkx2 matrix where m is number of LptD molecules, k is
%       number of vertices for each polygon and each entry is a row vector
%       [x,y] with position of each vertex

% OUTPUTS:
% NO argument outputs
% saves the model variable to a file called 'results.mat'

function model = mainModel(plotYes,settings,initPositions)

% we can input either settings or both settings and initial positions

if nargin == 2
    % we only have settings, not positions, so add the positions
    
    error('line 50 this code is not written')
    
elseif nargin < 2
    % first argument is settings always, this is all the data about
    % insertion rates and stuff
    
    settings.polygonSides = 251;

    settings.membraneCircumference = 500*pi; % um
    settings.currentMaxLen = 1000; %um
    settings.initialArea = settings.membraneCircumference*settings.currentMaxLen*2; % um^2

    settings.insRateProtein = 100; %um^2/s; estimate pi*0.0024*0.0024
    settings.insRateLPS = 100; % um^2/s

    settings.insRateBAM = 0; %s^-1; estimate 7/60
    settings.insRateLptD = 0; % per BAM

    settings.growthRate = settings.insRateLPS*settings.insRateLptD*settings.insRateBAM + settings.insRateProtein*settings.insRateBAM; % units I think um^2/s^2
    settings.sqrtGrowthRate = sqrt(settings.growthRate); % units I think um/s

    settings.BAMsize = 1;

    settings.time = 0;
    settings.dt = 0.01; %s
    settings.maxTime = 0.03;
    
    settings.surfaceTensionFlag = 1;
    settings.surfaceTensionStrength = 1;
    
    settings.splitFlag = 1;

    settings.proteinAddedNewInsertion = settings.insRateProtein*settings.dt;
    settings.LPSAddedNewInsertion = settings.insRateLPS*settings.dt;
    
    settings.saveLocation = '';
    settings.plotEvery = 25;
end

if nargin < 2
    % we have the settings, but not the initial positions of BAM and LptD.
    % Set these positions now
    
    initPositions.BAMlocs = [0,0.5*pi];
    initPositions.LptDlocs = [];
    initPositions.proteinVertices = {};
    initPositions.lpsVertices = {};
    %initPositions.lpsVertices = NaN(settings.polygonSides,2);
    
    % for all existing insertions add the corresponding material (note this
    % assumes one time step, if you want to skip you must pass the
    % vertices as an input)
    
    for i=1:size(initPositions.BAMlocs,1)
        initPositions.BAMlocs(i,:)
        vertices = findVerticesNewMaterialCircle(initPositions.BAMlocs(i,:),settings.polygonSides,0,settings.proteinAddedNewInsertion);
        initPositions.proteinVertices{i} = vertices;
        
        % note here we need to use settings.xyz as the plain non
        % settings.xyz parameters are set after this if statement
    end

    for i=1:size(initPositions.LptDlocs,1)
        vertices = findVerticesNewMaterialCircle(initPositions.LptDlocs(i,:),settings.polygonSides,0,settings.LPSAddedNewInsertion);
        initPositions.lpsVertices{i} = vertices;
    end
    
end

% set all the used parameters now, either to the values that were passed
% (if nargin >= 1) or to the values set above (if nargin < 1)
% we have to set the variables that are calculated from settings ones here
% in case we passed in settings

polygonSides = settings.polygonSides;

membraneCircumference = settings.membraneCircumference; % um
smVecs = [zeros(21,1),((-10:10)')*membraneCircumference]; % used in calcFlow.m

currentMaxLen = settings.currentMaxLen; %um
initialArea = settings.membraneCircumference*settings.currentMaxLen*2; % um^2
settings.initialArea = initialArea;

insRateProtein = settings.insRateProtein; %um^2/s; estimate pi*0.0024*0.0024
insRateLPS = settings.insRateLPS; % per LptD

insRateBAM = settings.insRateBAM; %s^-1; estimate 7/60
insRateLptD = settings.insRateLptD; % per BAM

growthRate = settings.insRateLPS*settings.insRateLptD*settings.insRateBAM + settings.insRateProtein*settings.insRateBAM; % units I think um^2/s^2
sqrtGrowthRate = sqrt(growthRate); % units I think um/s
settings.growthRate = growthRate;
settings.sqrtGrowthrate = sqrtGrowthRate;

BAMsize = settings.BAMsize;

surfaceTensionFlag = settings.surfaceTensionFlag;
surfaceTensionStrength = settings.surfaceTensionStrength;
splitFlag = settings.splitFlag;

time = settings.time;
dt = settings.dt; %s
maxTime = settings.maxTime;

proteinAddedNewInsertion = settings.insRateProtein*settings.dt;
settings.proteinAddedNewInsertion = proteinAddedNewInsertion;

LPSAddedNewInsertion = settings.insRateLPS*settings.dt;
settings.LPSAddedNewInsertion = LPSAddedNewInsertion;

saveLocation = settings.saveLocation;
plotEvery = settings.plotEvery;

% all information will be saved in the model struct

BAMlocs = initPositions.BAMlocs;
%BAMlocsInit = initPositions.BAMlocs; % store initial BAMlocs
LptDlocs = initPositions.LptDlocs;
%LptDlocsInit = initPositions.LptDlocs; % store initial LptDlocs
proteinVertices = initPositions.proteinVertices;
lpsVertices = initPositions.lpsVertices;

% save initial positions of everything
model.proteinVertices = proteinVertices;
model.lpsVertices = lpsVertices;
model.BAMlocs = BAMlocs;
model.LptDlocs = LptDlocs;

% store indices of insertion points which are actually relevant
proteinVerticesBAMs = {}; % This stores the indices of the BAMs that can actually affect each polygon
proteinVerticesLptDs = {};
lpsVerticesBAMs = {};
lpsVerticesLptDs = {};
%rightEdgeBAMs = [];
%leftEdgeBAMs = [];

% has to be a cell as can have different lengths. Guess it could have zeros
% but would need to think about this.

halfLen = currentMaxLen/2;
rightEdge = [halfLen*ones(100,1),linspace(0,membraneCircumference,100)'];
leftEdge = [-halfLen*ones(100,1),linspace(0,membraneCircumference,100)'];

model.rightEdge = rightEdge;
model.leftEdge = leftEdge;

% save the settings used to construct the model in model as well

model.settings = settings;

fig = figure;
visualiseSimple(BAMlocs,proteinVertices,LptDlocs,lpsVertices,rightEdge,leftEdge);
title('initial')
saveas(fig,[saveLocation,'initial.png']);
close(fig)

% insRateProtein
% insRateLPS
% proteinAddedNewInsertion
% LPSAddedNewInsertion

insertionAccuracy = 3;
maxInsDistProtein = (insRateProtein*exp(insertionAccuracy))/(2*pi);
maxInsDistLPS = (insRateProtein*exp(insertionAccuracy))/(2*pi);

maxInsDistProtein2 = maxInsDistProtein^2;
maxInsDistLPS2 = maxInsDistLPS^2;

% while time is less than max time
itCount = 0;
prematureEnd = 0;

% things we might want to calculate

cellLength = zeros(1,floor(maxTime/dt));
numBAMs = zeros(1,floor(maxTime/dt));
numLptDs = zeros(1,floor(maxTime/dt));

while time < maxTime && prematureEnd == 0
    
    itCount = itCount + 1;
    
    % calculate stuffs
    numBAMs(itCount) = size(BAMlocs,1);
    numLptDs(itCount) = size(LptDlocs,1);
    
    % calculate the mean left and right edge position
    
    meanLeftEdge = mean(leftEdge(:,1));
    meanRightEdge = mean(rightEdge(:,1));
    cellLength(itCount) = meanRightEdge - meanLeftEdge;
    
    % we need to maintain a list of newly inserted insertion points, or at
    % least their indices
    % there can only be at most one new bam location, and this is already
    % noted
    
    newBAMlocs = [];
    newLptDlocs = [];
    
    % add new LptD first so that new BAMs cannot insert new LptD
    
    % see if we need to add new lptD 
    
    newLptDIndex = 1;
    
    for BAM = 1:size(BAMlocs,1)
        % for each BAM we roll a die and if it is less than some number we
        % add an LptD near that BAM
        
        a = rand(1);
        
        %if count == 4
        if a < insRateLptD*dt
            % we pick a random angle in [0,2pi)
            theta = rand(1)*2*pi;
            
            % add twice BAMsize for space for LptD as well
            newLptDlocs(newLptDIndex,:) = [BAMlocs(BAM,1)+BAMsize*cos(theta),BAMlocs(BAM,2)+BAMsize*sin(theta)];
            
            newLptDIndex = newLptDIndex + 1;
        end
        
    end
    
    % add the new LptD to the LptD location list
    
    for lptD = 1:(newLptDIndex-1)
        newLptDLocsIndex = size(LptDlocs,1) + 1;

        LptDlocs(newLptDLocsIndex,:) = newLptDlocs(lptD,:);
        
        % check if this lptd should be considered for each protein polygon
        % add the new insertion point to all of the index lists for
        % each polygon
        for poly = 1:numel(proteinVerticesLptDs)
            % but only if the insertion point is close to at least one
            % vertex
            pos = proteinVertices{poly}(1,:);
            if sum((pos-newLptDlocs(lptD,:)).^2) < maxInsDistLPS2
                proteinVerticesLptDs{poly}(numel(proteinVerticesLptDs{poly})+1) = newLptDLocsIndex;
            end
        end
        for poly = 1:numel(lpsVerticesLptDs)
            % but only if the insertion point is close to at least one
            % vertex
            pos = lpsVertices{poly}(1,:);
            if sum((pos-newLptDlocs(lptD,:)).^2) < maxInsDistLPS2
                lpsVerticesLptDs{poly}(numel(lpsVerticesLptDs{poly})+1) = newLptDLocsIndex;
            end
        end
    end
    
    % see if we need to insert a new BAM
    
    newBamRandom = rand(1);
    
    % work out current area
    
    %currentArea = initialArea*exp(sqrtGrowthRate*time);
    
    if newBamRandom < insRateBAM*dt
    %if newBamRandom < insRateBAM*dt*currentArea
        % if yes add a BAM to a new randomly chosen location
        
        newBAMloc = rand(1,2); % gives two uniform random numbers
        
        % adjust the x one to be uniform between -currentMaxLen and
        % currentMaxLen
        % and y one to be uniform between 0 and membraneCircumference
        
        if settings.BAMtype == 0
            newBAMloc(1) = newBAMloc(1)*cellLength(itCount)-meanRightEdge;
            newBAMloc(2) = newBAMloc(2)*membraneCircumference;
        else
            newBAMloc(1) = newBAMloc(1)*500-250;
            newBAMloc(2) = newBAMloc(2)*membraneCircumference;
        end
        
        % we do not allow the insertion if it is within distance 1 of an
        % already existing BAM
        
        allowedInsertion = 1;
        
        for i=1:size(BAMlocs,1)
            if findDist(BAMlocs(i,:),[newBAMloc(1),newBAMloc(2)]) < model.settings.BAMsize
                allowedInsertion = 0;
            end
        end

        if allowedInsertion == 1
            newBAMIndex = size(BAMlocs,1) + 1;

            BAMlocs(newBAMIndex,:) = newBAMloc;

            newBAMlocs = newBAMloc;
            
            % add the new insertion point to all of the index lists for
            % each polygon
            for poly = 1:numel(proteinVerticesBAMs)
                % but only if the insertion point is close to at least one
                % vertex
                pos = proteinVertices{poly}(1,:);
                if sum((pos-newBAMloc).^2) < maxInsDistProtein2
                    proteinVerticesBAMs{poly}(numel(proteinVerticesBAMs{poly})+1) = newBAMIndex;
                end
            end
            for poly = 1:numel(lpsVerticesBAMs)
                % but only if the insertion point is close to at least one
                % vertex
                pos = lpsVertices{poly}(1,:);
                if sum((pos-newBAMloc).^2) < maxInsDistProtein2
                    lpsVerticesBAMs{poly}(numel(lpsVerticesBAMs{poly})+1) = newBAMIndex;
                end
            end
            
        end
        
    end
    
    % move protein polygon vertices
    % have to loop through all polygons, then all pairs of vertices
    for poly = 1:numel(proteinVertices)
        if numel(LptDlocs) > 0
            flow = calcFlow(proteinVertices{poly},BAMlocs(proteinVerticesBAMs{poly},:),insRateProtein,LptDlocs(proteinVerticesLptDs{poly},:),insRateLPS,smVecs);
        else
            flow = calcFlow(proteinVertices{poly},BAMlocs(proteinVerticesBAMs{poly},:),insRateProtein,[],insRateLPS,smVecs);
        end
        proteinVertices{poly} = proteinVertices{poly} + flow*dt;
    end

    % move lps polygon vertices
    % have to loop through all polygons, then all pairs of vertices

    for poly = 1:numel(lpsVertices)
        flow = calcFlow(lpsVertices{poly},BAMlocs(lpsVerticesBAMs{poly},:),insRateProtein,LptDlocs(lpsVerticesLptDs{poly},:),insRateLPS,smVecs);
        lpsVertices{poly} = lpsVertices{poly} + flow*dt;
    end
    
    % move the edges
    
    flows = calcFlow(rightEdge,BAMlocs,insRateProtein,LptDlocs,insRateLPS,smVecs);
    rightEdge = rightEdge + flows*dt;
    
    flows = calcFlow(leftEdge,BAMlocs,insRateProtein,LptDlocs,insRateLPS,smVecs);
    leftEdge = leftEdge + flows*dt;

    % move insertion points (lps)
    
    % have to record new positions separately and move everything
    % simultaneously
    movedLptDlocs = zeros(size(LptDlocs));
    
    for j = 1:size(LptDlocs,1)
        % find the flow from all points EXCEPT itself
        tempLptDlocs = LptDlocs;
        tempLptDlocs(j,:) = [];
        
        flow = calcFlow(LptDlocs(j,:),BAMlocs,insRateProtein,tempLptDlocs,insRateLPS,smVecs);

        movedLptDlocs(j,:) = LptDlocs(j,:) + flow*dt;
        
    end

    % move insertion points (protein)
    
    % we need to record the new positions separately and then move
    % everything simultaneously
    movedBAMlocs = zeros(size(BAMlocs));
    
    % this has to be done individually as different sets of BAMs are
    % allowed for each (cant move due to flow from itself)
    for j = 1:size(BAMlocs,1)
        % find the flow from all points EXCEPT itself
        tempBAMlocs = BAMlocs;
        tempBAMlocs(j,:) = [];

        flow = calcFlow(BAMlocs(j,:),tempBAMlocs,insRateProtein,LptDlocs,insRateLPS,smVecs);
        
        movedBAMlocs(j,:) = BAMlocs(j,:) + flow*dt;        
    end

    % now that we have the new locations for both BAM and LptD, record them
    
    LptDlocs = movedLptDlocs;
    BAMlocs = movedBAMlocs;

    % add new material from insertions
    
    % add new material (protein)
    % check to see whether new BAMs were added or not
    if ~isempty(newBAMlocs)
        % loop through them (in case somehow we have more than one added)
        for j=1:size(newBAMlocs,1)
            newIndex = numel(proteinVertices)+1;
            vertices = findVerticesNewMaterialCircle(BAMlocs(end-j+1,:),polygonSides,0,proteinAddedNewInsertion);
            proteinVertices{newIndex}(:,:) = vertices;
            
            % note down the insertion points that can actually affect this
            % polygon
            insertionPoints = [];
            index = 1;
            for k=1:size(BAMlocs,1)
                dists = sum((vertices - BAMlocs(k,:)).^2,2);
                
                % 22500 = 150^2
                closeDists = sum(dists < maxInsDistProtein2);
                if closeDists > 0
                    % store indices
                    insertionPoints(index) = k;
                    index = index + 1;
                end
            end
            
            proteinVerticesBAMs{newIndex} = insertionPoints;
            
            % do same for LPS
            insertionPoints = [];
            index = 1;
            for k=1:size(LptDlocs,1)
                dists = sum((vertices - LptDlocs(k,:)).^2,2);
                
                % 22500 = 150^2
                closeDists = sum(dists < maxInsDistLPS2);
                if closeDists > 0
                    % store indices
                    insertionPoints(index) = k;
                    index = index + 1;
                end
            end
            proteinVerticesLptDs{newIndex} = insertionPoints;
        end
    end
    
    % add new material (LPS)
    % check to see whether new LptDs were added or not
    if ~isempty(newLptDlocs)
        % loop through them (in case somehow we have more than one added)
        for j=1:size(newLptDlocs,1)
            newIndex = numel(lpsVertices)+1;
            
            % we want to use the moved position of the new lptd, not the
            % original position, so use the LptDlocs(end-j+1)
            vertices = findVerticesNewMaterialCircle(LptDlocs(end-j+1,:),polygonSides,0,LPSAddedNewInsertion);

            lpsVertices{newIndex}(:,:) = vertices;
            
            % note down the insertion points that can actually affect this
            % polygon
            insertionPoints = [];
            index = 1;
            for k=1:size(BAMlocs,1)
                dists = sum((vertices - BAMlocs(k,:)).^2,2);
                
                % 22500 = 150^2
                closeDists = sum(dists < maxInsDistProtein2);
                if closeDists > 0
                    % store indices
                    insertionPoints(index) = k;
                    index = index + 1;
                end
            end
            
            lpsVerticesBAMs{newIndex} = insertionPoints;
            
            % do same for LPS
            insertionPoints = [];
            index = 1;
            for k=1:size(LptDlocs,1)
                dists = sum((vertices - LptDlocs(k,:)).^2,2);
                
                % 22500 = 150^2
                closeDists = sum(dists < maxInsDistLPS2);
                if closeDists > 0
                    % store indices
                    insertionPoints(index) = k;
                    index = index + 1;
                end
            end
            lpsVerticesLptDs{newIndex} = insertionPoints;
        end
    end
    
        
    
    % move points under surface tension if required
    
    if surfaceTensionFlag == 1
        
        % we want to move the vertices of each polygon due to surface
        % tension
        
        % do protein first
        
        for poly = 1:numel(proteinVertices)
            
            for i=1:5
                newPoints = surfaceTensionPolygon(proteinVertices{poly}(:,:),dt,surfaceTensionStrength);
            end
            
            proteinVertices{poly}(:,:) = newPoints;
        end
        
        % then LPS, but only if lps vertices exist
        
        if numel(LptDlocs) > 0
            for poly = 1:numel(lpsVertices)
                
                for i=1:5
                    newPoints = surfaceTensionPolygon(lpsVertices{poly}(:,:),dt,surfaceTensionStrength);
                end
                
                lpsVertices{poly}(:,:) = newPoints;
            end
        end
        
    end
    
    % check for vertices being 'too close'
    % note that the variable 'tooClose', here set to 1, should be replaced
    % at some point to become a model parameter
    
    plotSplit = 0;
    if splitFlag == 1
        
        if mod(itCount,10) == 0
    
            [problems,newProteinVertices,indexRemoved] = checkPolygonDistances2(proteinVertices,1,0);

            if problems == 1
                disp('some regions too close, changing model')
                disp(['iteration: ',num2str(itCount)])

                % we have regions that are too close
                % snap a picture
                if plotSplit == 1
                    fig = figure;
                    visualiseSimple(model)

                    title(['some regions are too close',num2str(time)])
                    saveas(fig,[saveLocation,'tooClose-it',num2str(itCount),'.png']);
                    close(fig)
                end

                % resolve the issues
               
                proteinVertices = newProteinVertices;
                proteinVerticesBAMs(indexRemoved) = [];
                proteinVerticesLptDs(indexRemoved) = [];

                for i=1:numel(proteinVertices)
                    assert(size(proteinVertices,1) ~= 0)
                end

                % take a new picture of the resolution

                if plotSplit == 1
                    fig = figure;
                    visualiseSimple(model)

                    title(['resolution',num2str(time)])
                    saveas(fig,[saveLocation,'resolution-it',num2str(itCount),'.png']);
                    close(fig)
                end
            end
        end
    end
    
    % every so many iterations re-check how close vertices are to
    % insertions
    
    if mod(itCount,25) == 0
        for poly = 1:numel(proteinVertices)
            % PROTEIN VS BAM
            vertices = proteinVertices{poly};
            % note down the insertion points that can actually affect this
            % polygon
            insertionPoints = [];
            index = 1;
            for k=1:size(BAMlocs,1)
                dists = sum((vertices - BAMlocs(k,:)).^2,2);

                % 22500 = 150^2
                closeDists = sum(dists < maxInsDistProtein2);
                if closeDists > 0
                    % store indices
                    insertionPoints(index) = k;
                    index = index + 1;
                end
            end

            proteinVerticesBAMs{poly} = insertionPoints;
            
            % PROTEIN VS LPTD
            insertionPoints = [];
            index = 1;
            for k=1:size(LptDlocs,1)
                dists = sum((vertices - LptDlocs(k,:)).^2,2);

                % 22500 = 150^2
                closeDists = sum(dists < maxInsDistLPS2);
                if closeDists > 0
                    % store indices
                    insertionPoints(index) = k;
                    index = index + 1;
                end
            end

            proteinVerticesLptDs{poly} = insertionPoints;
        end
        % check lps too
        for poly = 1:numel(lpsVertices)
            % LPS VS BAM
            vertices = lpsVertices{poly};
            % note down the insertion points that can actually affect this
            % polygon
            insertionPoints = [];
            index = 1;
            for k=1:size(BAMlocs,1)
                dists = sum((vertices - BAMlocs(k,:)).^2,2);

                % 22500 = 150^2
                closeDists = sum(dists < maxInsDistProtein2);
                if closeDists > 0
                    % store indices
                    insertionPoints(index) = k;
                    index = index + 1;
                end
            end

            lpsVerticesBAMs{poly} = insertionPoints;
            
            % LPS VS LPTD
            vertices = lpsVertices{poly};
            % note down the insertion points that can actually affect this
            % polygon
            insertionPoints = [];
            index = 1;
            for k=1:size(LptDlocs,1)
                dists = sum((vertices - LptDlocs(k,:)).^2,2);

                % 22500 = 150^2
                closeDists = sum(dists < maxInsDistLPS2);
                if closeDists > 0
                    % store indices
                    insertionPoints(index) = k;
                    index = index + 1;
                end
            end

            lpsVerticesLptDs{poly} = insertionPoints;
           
        end
    end
    
    % take a picture
    if mod(itCount,plotEvery) == 0
        fig = figure;
        visualiseSimple(BAMlocs,proteinVertices,LptDlocs,lpsVertices,rightEdge,leftEdge);
        title(['time is ',num2str(time)])
        saveas(fig,[saveLocation,'it-',num2str(itCount),'.png']);
        %close(fig)
    end
    
    
    % calculate anything you want to calculate
    
    
    % SHRINK THE MEMBRANE
    
    % check to see if we have exceeded twice the original length
    
        % if we have, cut the left side out
        
        % first take a picture of how it looked
        
        % all polygons/BAMS/LptD with all vertices greater than zero are removed
        
        % all polygons/BAMS/LptD with a vertex greater than zero have that vertex set
        % to zero
        
        % all polygons/BAMS/LptD with a vertex less than 2L have that vertex
        % set to 2L
        
        % move all remaining vertices/BAMS/LptD to the right by L
        
        % take a picture of how it looked
    
    
    % increase time
    
    time = time + dt;
    
    % plot
    
    %visualiseSimple(model)

end

% save final positions of everything
model.proteinVertices = proteinVertices;
model.lpsVertices = lpsVertices;
model.BAMlocs = BAMlocs;
model.LptDlocs = LptDlocs;

model.rightEdge = rightEdge;
model.leftEdge = leftEdge;

model.cellLength = cellLength;
model.numBAMs = numBAMs;
model.numLptDs = numLptDs;

if plotYes == 1
    %visualiseSimple(model)

    fig = figure;
    visualiseSimple(BAMlocs,proteinVertices,LptDlocs,lpsVertices,rightEdge,leftEdge)
    
    % display a circle with same area as the protein region in black
    
    areaProtein = polyarea(proteinVertices{1}(:,1),proteinVertices{1}(:,2));
    radius = sqrt(areaProtein/pi);
    theta = linspace(0,2*pi,1000);
    x = BAMlocs(1,1) + radius*cos(theta);
    y = BAMlocs(1,2) + radius*sin(theta);
    
    plot(x,y,'k--')

    title(['reached end, time is ',num2str(time)])
    saveas(fig,[saveLocation,'endPoint.png']);

    close(fig)

end

save([saveLocation,'results.mat'],'model')

end