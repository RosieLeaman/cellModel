function model = runModel(m,useGivenInsertions)
% m is number of iterations
% useGivenInsertions is 0 if not and a file location if yes

if nargin < 1
    m = 3000;
end

if useGivenInsertions ~= 0
    data = load(useGivenInsertions);
    oldModel = data.model;
    BAMlocs = [];
    index = 1;
    for j=1:size(oldModel.BAMlocs,1)
        if oldModel.BAMlocs(j,1) > -2000 && oldModel.BAMlocs(j,1) < 0
            BAMlocs(index,:) = oldModel.BAMlocs(j,:);
            index = index + 1;
        end
    end
    BAMlocs(:,1) = BAMlocs(:,1) + 1000;

    LptDlocs = [];
    index = 1;
    for j=1:size(oldModel.LptDlocs,1)
        if oldModel.LptDlocs(j,1) > -2000 && oldModel.LptDlocs(j,1) < 0
            LptDlocs(index,:) = oldModel.LptDlocs(j,:);
            index = index + 1;
        end
    end
    LptDlocs(:,1) = LptDlocs(:,1) + 1000;
    
    initPositions.BAMlocs = BAMlocs;
    
    initPositions.LptDlocs = LptDlocs;
else
    initPositions.BAMlocs = [];
    
    initPositions.LptDlocs = [];
end

% set the random seed and record it
rng('shuffle'); % first re-set the random generator, otherwise get same results each time
randSeed = rand(1); % pick a random seed
rng(randSeed); % re-seed with this seed
settings.randSeed = randSeed;

% set up non-changing settings

settings.polygonSides = 250;
settings.membraneCircumference = 500*pi;
settings.currentMaxLen = 2000;
settings.BAMsize = 1;
settings.time = 0;
settings.dt = 0.01;
settings.maxTime = m*0.01;
settings.surfaceTensionFlag = 0;
settings.surfaceTensionStrength = 0;
settings.splitFlag = 1;
settings.plotEvery = 100;

settings.insRateProtein = 1500;
settings.insRateLPS = 1500;
settings.insRateBAM = 10;
settings.insRateLptD = 0.1;  % 0.022 = 300/15(30^3) p46 lab book 4

settings.BAMtype = 0; %0 = random 1 = mid-cell

% set up save location and file structures
% find our current directory
currentFolder = pwd();

% make a new directory with current date and time
time = clock; % [year - month - day - hour - min - s]
timeString = [num2str(time(1)),'-',num2str(time(2)),'-',num2str(time(3))];
timeString = [timeString,'-',num2str(time(4)),'-',num2str(time(5))];
mkdir([currentFolder,'/test',timeString]);
bigSaveFolder = [currentFolder,'/test',timeString];

% set up initial BAM locations
% set up initial LptD inside loop as we want to vary this later

% If there are initial locations we only send the location of BAMs/LptDs so
% that the correct insertion rate is used for the initial protein vertices
% see if they were passed 

for i = [1500,3000]
    i
    % set up the settings that vary
    settings.insRateProtein = i;
    
    settings.saveLocation = [bigSaveFolder,'/test-proteinInsertionRate',num2str(i)];

    % run the model with given settings
    
    %try
    model = mainModel(0,settings,initPositions);
    
    figure;
    visualiseSimple(model.BAMlocs,model.proteinVertices,model.LptDlocs,model.lpsVertices,model.rightEdge,model.leftEdge);
    %catch e
%         e
%         visualiseSimple(model.BAMlocs,model.proteinVertices,model.LptDlocs,model.lpsVertices,model.rightEdge,model.leftEdge);
%     end
end
