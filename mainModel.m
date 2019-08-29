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
    
    initPositions.BAMlocs = [[0,0.5*pi]];
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

model.BAMlocs = initPositions.BAMlocs;
model.BAMlocsInit = initPositions.BAMlocs; % store initial BAMlocs
model.LptDlocs = initPositions.LptDlocs;
model.LptDlocsInit = initPositions.LptDlocs; % store initial LptDlocs
model.proteinVertices = initPositions.proteinVertices;
model.lpsVertices = initPositions.lpsVertices;

% save the settings used to construct the model in model as well

model.settings = settings;

% while time is less than max time
count = 0;
prematureEnd = 0;

testFlag = 0;

insRateProtein
insRateLPS
proteinAddedNewInsertion
LPSAddedNewInsertion

while time < maxTime && prematureEnd == 0
    
    count = count + 1;
    
    if testFlag == 1
        disp(['COUNT IS: ',count])
        model.proteinVertices
    end

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
        % if yes add a BAM to a new randomly chosen location
        
        newBAMloc = rand(1,2); % gives two uniform random numbers
        
        % adjust the x one to be uniform between -currentMaxLen and
        % currentMaxLen
        % and y one to be uniform between 0 and membraneCircumference
        newBAMloc(1) = newBAMloc(1)*2*currentMaxLen-currentMaxLen;
        newBAMloc(2) = newBAMloc(2)*membraneCircumference;
        
        % we do not allow the insertion if it is within distance 1 of an
        % already existing BAM
        
        allowedInsertion = 1;
        
        for i=1:size(model.BAMlocs,1)
            if findDist(model.BAMlocs(i,:),[newBAMloc(1),newBAMloc(2)]) < model.settings.BAMsize
                allowedInsertion = 0;
            end
        end
        
        if allowedInsertion == 1
            newBAMIndex = size(model.BAMlocs,1) + 1;

            model.BAMlocs(newBAMIndex,:) = newBAMloc;

            newBAMlocs = newBAMloc;
        end
        
    end
       
    % see if we need to add new lptD 
    
    newLptDIndex = 1;
    
    for BAM = 1:size(model.BAMlocs,1)
        % for each BAM we roll a die and if it is less than some number we
        % add an LptD near that BAM
        
        a = rand(1);
        
        %if count == 4
        if a < insRateLptD*dt
            % we pick a random angle in [0,2pi)
            theta = rand(1)*2*pi;
            
            % add twice BAMsize for space for LptD as well
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
    for poly = 1:numel(model.proteinVertices)
        for j=1:size(model.proteinVertices{poly},1)

            flow = calcFlow(model.proteinVertices{poly}(j,:),model.BAMlocs,insRateProtein,model.LptDlocs,insRateLPS,membraneCircumference,0);
            
            newPos = findNewVertexPosition(model.proteinVertices{poly}(j,:),flow,dt,membraneCircumference);
            
            model.proteinVertices{poly}(j,:) = newPos;
        end
    end
    
    if testFlag == 1
        disp('AFTER MOVING POLYGONS IS: ')
        model.proteinVertices
    end
    
    % move lps polygon vertices
    % have to loop through all polygons, then all pairs of vertices

    for poly = 1:numel(model.lpsVertices)
        for j=1:size(model.lpsVertices{poly},1)

            flow = calcFlow(model.lpsVertices{poly}(j,:),model.BAMlocs,insRateProtein,model.LptDlocs,insRateLPS,membraneCircumference,0);
            
            newPos = findNewVertexPosition(model.lpsVertices{poly}(j,:),flow,dt,membraneCircumference);
            
            model.lpsVertices{poly}(j,:) = newPos;
        end
    end
    
    if testFlag == 1
        disp('AFTER MOVING LPS IS: ')
        model.proteinVertices
    end
    
    
    % move insertion points (lps)
    
    % have to record new positions separately and move everything
    % simultaneously
    movedLptDlocs = zeros(size(model.LptDlocs));
    
    for j = 1:size(model.LptDlocs,1)
        % find the flow from all points EXCEPT itself
        tempLptDlocs = model.LptDlocs;
        tempLptDlocs(j,:) = [];
        
        flow = calcFlow(model.LptDlocs(j,:),model.BAMlocs,insRateProtein,tempLptDlocs,insRateLPS,membraneCircumference,0);
            
        newPos = findNewVertexPosition(model.LptDlocs(j,:),flow,dt,membraneCircumference);
            
        %model.LptDlocs(j,:) = newPos;
        movedLptDlocs(j,:) = newPos;
        
    end
    
    if testFlag == 1
        disp('AFTER MOVING LPTD IS: ')
        model.proteinVertices
    end
    
    % move insertion points (protein)
    
    % we need to record the new positions separately and then move
    % everything simultaneously
    movedBAMlocs = zeros(size(model.BAMlocs));
        
    for j = 1:size(model.BAMlocs,1)
        % find the flow from all points EXCEPT itself
        tempBAMlocs = model.BAMlocs;
        tempBAMlocs(j,:) = [];

        flow = calcFlow(model.BAMlocs(j,:),tempBAMlocs,insRateProtein,model.LptDlocs,insRateLPS,membraneCircumference,0);
        
        
        newPos = findNewVertexPosition(model.BAMlocs(j,:),flow,dt,membraneCircumference);
            
        %model.BAMlocs(j,:) = newPos;
        movedBAMlocs(j,:) = newPos;
        
    end
    
    if testFlag == 1
        disp('AFTER MOVING BAMS IS: ')
        model.proteinVertices
    end
    
    % now that we have the new locations for both BAM and LptD, record them
    
    model.LptDlocs = movedLptDlocs;
    model.BAMlocs = movedBAMlocs;

    % add new material from insertions
    
    % add new material (protein)
    % check to see whether new BAMs were added or not
    if ~isempty(newBAMlocs)
        % loop through them (in case somehow we have more than one added)
        for j=1:size(newBAMlocs,1)
            newIndex = numel(model.proteinVertices)+1;
            vertices = findVerticesNewMaterialCircle(model.BAMlocs(end-j+1,:),polygonSides,0,proteinAddedNewInsertion);
            model.proteinVertices{newIndex}(:,:) = vertices;

        end
    end
    
    if testFlag == 1
        disp('AFTER ADDING NEW BAMS IS: ')
        model.proteinVertices
    end
    
    % add new material (LPS)
    % check to see whether new LptDs were added or not
    if ~isempty(newLptDlocs)
        % loop through them (in case somehow we have more than one added)
        for j=1:size(newLptDlocs,1)
            newIndex = numel(model.lpsVertices)+1;
            
            % we want to use the moved position of the new lptd, not the
            % original position, so use the model.LptDlocs(end-j+1)
            vertices = findVerticesNewMaterialCircle(model.LptDlocs(end-j+1,:),polygonSides,0,LPSAddedNewInsertion);

            model.lpsVertices{newIndex}(:,:) = vertices;
        end
    end
        
    
    % move points under surface tension if required
    
    if settings.surfaceTensionFlag == 1
        
        % we want to move the vertices of each polygon due to surface
        % tension
        
        % do protein first
        
        for poly = 1:numel(model.proteinVertices)
            
            for i=1:5
                newPoints = surfaceTensionPolygon(model.proteinVertices{poly}(:,:),dt,surfaceTensionStrength);
            end
            
            model.proteinVertices{poly}(:,:) = newPoints;
        end
        
        % then LPS, but only if lps vertices exist
        
        if numel(model.LptDlocs) > 0
            for poly = 1:numel(model.lpsVertices)
                
                for i=1:5
                    newPoints = surfaceTensionPolygon(model.lpsVertices{poly}(:,:),dt,surfaceTensionStrength);
                end
                
                model.lpsVertices{poly}(:,:) = newPoints;
            end
        end
        
    end
    
    % check for vertices being 'too close'
    % note that the variable 'tooClose', here set to 1, should be replaced
    % at some point to become a model parameter
    
    plotSplit = 0;
    if splitFlag == 1
    
        [problems,newModel] = checkPolygonDistances2(model,1,0);

        if problems == 1
            %disp('we have located a problem')
            % we have regions that are too close
            % snap a picture
            if plotSplit == 1
                fig = figure;
                visualiseSimple(model)

                title(['some regions are too close',num2str(time)])
                saveas(fig,[saveLocation,'tooClose-it',num2str(count),'.png']);
                close(fig)
            end
        
            % resolve the issues

            model = newModel;
        
            % take a new picture of the resolution

            if plotSplit == 1
                fig = figure;
                visualiseSimple(model)

                title(['resolution',num2str(time)])
                saveas(fig,[saveLocation,'resolution-it',num2str(count),'.png']);
                close(fig)
            end

        end
    end
    
    % take a picture
    if mod(count,plotEvery) == 0
        fig = figure;
        visualiseSimple(model);
        title(['time is ',num2str(time)])
        saveas(fig,[saveLocation,'it-',num2str(count),'.png']);
        close(fig)
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

if plotYes == 1
    %visualiseSimple(model)

    fig = figure;
    visualiseSimple(model)
    
    % display a circle with same area as the protein region in black
    
    areaProtein = polyarea(model.proteinVertices{1}(:,1),model.proteinVertices{1}(:,2));
    radius = sqrt(areaProtein/pi);
    theta = linspace(0,2*pi,1000);
    x = model.BAMlocs(1,1) + radius*cos(theta);
    y = model.BAMlocs(1,2) + radius*sin(theta);
    
    plot(x,y,'k--')

    title(['reached end, time is ',num2str(time)])
    saveas(fig,[saveLocation,'endPoint.png']);

    close(fig)

end


save([saveLocation,'results.mat'],'model')