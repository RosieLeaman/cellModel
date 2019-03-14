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


function mainModel(plotYes,settings,initPositions)

% we can input either settings or both settings and initial positions

if nargin < 2
    % first argument is settings always, this is all the data about
    % insertion rates and stuff
    
    settings.polygonSides = 251;

    settings.membraneCircumference = pi; % um
    settings.currentMaxLen = 2; %um
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
    settings.maxTime = 0.1;
    
    settings.surfaceTensionFlag = 1;

    settings.proteinAddedNewInsertion = settings.insRateProtein*settings.dt;
    settings.LPSAddedNewInsertion = settings.insRateLPS*settings.dt;
    
    settings.saveLocation = '';

elseif nargin < 2
    % we have the settings, but not the initial positions of BAM and LptD.
    % Set these positions now
    
    initPositions.BAMlocs = [[0,0.5*pi]];
    initPositions.LptDlocs = [];
    initPositions.proteinVertices = [];
    initPositions.lpsVertices = NaN(settings.polygonSides,2);
    
    % for all existing insertions add the corresponding material (note this
    % assumes one time step, if you want to skip you must pass the
    % vertices as an input)
    
    for i=1:size(initPositions.BAMlocs,1)
        initPositions.BAMlocs(i,:)
        vertices = findVerticesNewMaterialCircle(initPositions.BAMlocs(i,:),settings.polygonSides,settings.proteinAddedNewInsertion);
        initPositions.proteinVertices(:,:,i) = vertices;
        
        % note here we need to use settings.xyz as the plain non
        % settings.xyz parameters are set after this if statement
    end

    for i=1:size(initPositions.LptDlocs,1)
        vertices = findVerticesNewMaterial(initPositions.LptDlocs(i,:),settings.polygonSides,settings.LPSAddedNewInsertion);
        initPositions.lpsVertices(:,:,i) = vertices;
    end
    
end

% set all the used parameters now, either to the values that were passed
% (if nargin >= 1) or to the values set above (if nargin < 1)

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

proteinAddedNewInsertion = settings.proteinAddedNewInsertion;
LPSAddedNewInsertion = settings.LPSAddedNewInsertion;

% all information will be saved in the model struct

model.BAMlocs = initPositions.BAMlocs;
model.BAMlocsInit = initPositions.BAMlocs; % store initial BAMlocs
model.LptDlocs = initPositions.LptDlocs;
model.LptDlocsInit = initPositions.LptDlocs; % store initial LptDlocs
model.proteinVertices = initPositions.proteinVertices;
model.lpsVertices = initPositions.lpsVertices;

% save the settings used to construct the model in model as well

model.settings = settings;

% make a pretty plot (if requested)
if plotYes == 1
    figure;hold on;

    visualiseSimple(model)
end


% while time is less than max time
count = 0;
while time < maxTime
    count = count + 1;

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
        
        %if count == 4
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

            flow = calcFlow(model.proteinVertices(j,:,poly),model.BAMlocs,insRateProtein,model.LptDlocs,insRateLPS,membraneCircumference,0);
            
            newPos = findNewVertexPosition(model.proteinVertices(j,:,poly),flow,dt,membraneCircumference);
            
            model.proteinVertices(j,:,poly) = newPos;
        end
    end
    
    % move lps polygon vertices
    % have to loop through all polygons, then all pairs of vertices

    for poly = 1:size(model.lpsVertices,3)
        for j=1:size(model.lpsVertices(:,:,1),1)

            flow = calcFlow(model.lpsVertices(j,:,poly),model.BAMlocs,insRateProtein,model.LptDlocs,insRateLPS,membraneCircumference,0);
            
            newPos = findNewVertexPosition(model.lpsVertices(j,:,poly),flow,dt,membraneCircumference);
            
            model.lpsVertices(j,:,poly) = newPos;
        end
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
    
    % now that we have the new locations for both BAM and LptD, record them
    
    model.LptDlocs = movedLptDlocs;
    model.BAMlocs = movedBAMlocs;

    % add new material from insertions
    
    % add new material (protein)
    % check to see whether new BAMs were added or not
    if ~isempty(newBAMlocs)
        % loop through them (in case somehow we have more than one added)
        for j=1:size(newBAMlocs,1)
            vertices = findVerticesNewMaterial(newBAMlocs(j,:),polygonSides,proteinAddedNewInsertion);
            model.proteinVertices(:,:,newBAMIndex) = vertices;

            if plotYes == 1
                plot(newBAMlocs(j,1),newBAMlocs(j,2),'xk','linewidth',2)
            end
           
        end
    end
    
    % add new material (LPS)
    % check to see whether new LptDs were added or not
    if ~isempty(newLptDlocs)
        % loop through them (in case somehow we have more than one added)
        for j=1:size(newLptDlocs,1)
            vertices = findVerticesNewMaterial(newLptDlocs(j,:),polygonSides,LPSAddedNewInsertion);
            
            % test whether we already have some LPS vertices or not
            if isnan(model.lpsVertices(1,1))
                newLptDIndex = 1;
            else
                newLptDIndex = size(model.lpsVertices,3) + 1;
            end
            
            model.lpsVertices(:,:,newLptDIndex) = vertices;

            if plotYes == 1
                plot(newLptDlocs(j,1),newLptDlocs(j,2),'dk','linewidth',2)
            end
           
        end
    end
        
    
    % move points under surface tension if required
    
    if settings.surfaceTensionFlag == 1
        
        % we want to move the vertices of each polygon due to surface
        % tension
        
        % do protein first
        
        for poly = 1:size(model.proteinVertices,3)
            
            newPoints = surfaceTensionPolygon(model.proteinVertices(:,:,poly),dt);
            
            model.proteinVertices(:,:,poly) = newPoints;
        end
        
        % then LPS, but only if lps vertices exist
        
        if ~isempty(model.LptDlocs)
            for poly = 1:size(model.lpsVertices,3)

                newPoints = surfaceTensionPolygon(model.lpsVertices(:,:,poly),dt);

                model.lpsVertices(:,:,poly) = newPoints;
            end
        end
        
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
    visualiseSimple(model)

    figure;
    visualiseSimple(model)

    visualise(model)
end


save([saveLocation,'results.mat'],'model')