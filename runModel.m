function runModel()

% set up non-changing settings

settings.polygonSides = 1000;
settings.membraneCircumference = 500*pi;
settings.currentMaxLen = 1000;
settings.BAMsize = 1;
settings.time = 0;
settings.dt = 0.01;
settings.maxTime = 5;
settings.surfaceTensionFlag = 0;
settings.surfaceTensionStrength = 0;
settings.splitFlag = 1;
settings.plotEvery = 25;

settings.insRateProtein = 100;
settings.insRateLPS = 100;
settings.insRateBAM = 0.1;
settings.insRateLptD = 0.1;

% test with different BAM insertion strengths

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

% set up initial locations

initPositions.BAMlocs = [[100,100]];

% the amount of material determines the size of the circle as amount =
% pi*r^2
radius = 2;
vertices = findVerticesNewMaterialCircle(initPositions.BAMlocs(1,:),settings.polygonSides,1,radius);
initPositions.proteinVertices{1} = vertices;

initPositions.LptDlocs = [];
initPositions.lpsVertices = {};

for i = 10000
    i
    % set up the settings that vary
    settings.insRateProtein = i;
    
    settings.saveLocation = [bigSaveFolder,'/test-proteinInsertionRate',num2str(i)];

    % run the model with given settings

    mainModel(1,settings,initPositions);
end
