function runModel()

% set up save location and file structures

% find our current directory
currentFolder = pwd();

% make a new directory with current date and time
time = clock; % [year - month - day - hour - min - s]
timeString = [num2str(time(1)),'-',num2str(time(2)),'-',num2str(time(3))];
timeString = [timeString,'-',num2str(time(4)),'-',num2str(time(5))];
mkdir([currentFolder,'/test',timeString]);
bigSaveFolder = [currentFolder,'/test',timeString];

% set up non-changing settings

settings.polygonSides = 251;
settings.membraneCircumference = 500*pi;
settings.currentMaxLen = 1000;
settings.BAMsize = 1;
settings.time = 0;
settings.dt = 0.01;
settings.maxTime = 0.01;
settings.surfaceTensionFlag = 1;

% set up initial BAM locations
% set up initial LptD inside loop as we want to vary this later

% set up initial locations

initPositions.BAMlocs = [[100,100]];

% the amount of material determines the size of the circle as amount =
% pi*r^2
radius = 2;
vertices = findVerticesNewMaterialCircle(initPositions.BAMlocs(1,:),settings.polygonSides,pi*radius^2);
initPositions.proteinVertices(:,:,1) = vertices;

for i = linspace(10,30,3)
    % set up the settings that vary
    settings.insRateProtein = i;
    settings.insRateLPS = 10;
    settings.insRateBAM = 0;
    settings.insRateLptD = 0;
  
    settings.saveLocation = [bigSaveFolder,'/test-proteinInsertionRate',num2str(i)];

    initPositions.LptDlocs = [];
    initPositions.lpsVertices = NaN(settings.polygonSides,2);

    % run the model with given settings

    mainModel(1,settings,initPositions);

end
