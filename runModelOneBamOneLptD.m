function runModelOneBamOneLptD()

% set up save location and file structures

% find our current directory
currentFolder = pwd();

% make a new directory with current date and time
time = clock; % [year - month - day - hour - min - s]
timeString = [num2str(time(1)),'-',num2str(time(2)),'-',num2str(time(3))];
timeString = [timeString,'-',num2str(time(4)),'-',num2str(time(5))];
mkdir([currentFolder,'/oneBamOneLptD',timeString]);
bigSaveFolder = [currentFolder,'/oneBamOneLptD',timeString];

% set up non-changing settings

settings.polygonSides = 251;
settings.membraneCircumference = 500*pi;
settings.currentMaxLen = 1000;
settings.BAMsize = 1;
settings.time = 0;
settings.dt = 0.01;
settings.maxTime = 2;
settings.surfaceTensionFlag = 1;

% set up initial BAM locations
% set up initial LptD inside loop as we want to vary this later

% set up initial locations

initPositions.BAMlocs = [[100,100]];

% the amount of material determines the size of the circle as amount =
% pi*r^2
radius = 6;
vertices = findVerticesNewMaterialCircle(initPositions.BAMlocs(1,:),settings.polygonSides,0,pi*radius^2);
initPositions.proteinVertices{1}(:,:) = vertices;

for i = linspace(5.1,5.1,1)
    % set up size of protein
    disp(i)
    
    radius = 12.1;
    vertices = findVerticesNewMaterialCircle(initPositions.BAMlocs(1,:),settings.polygonSides,0,pi*radius^2);
    initPositions.proteinVertices{1}(:,:) = vertices;
    
    % set up the settings that vary
    settings.insRateProtein = 250;
    settings.insRateLPS = 250;
    settings.insRateBAM = 0;
    settings.insRateLptD = 0;
    settings.surfaceTensionStrength = 100;
    settings.surfaceTensionFlag = 1;
  
    settings.saveLocation = [bigSaveFolder,'/test-proteinInsertionRate',num2str(i)];

    initPositions.LptDlocs = [108,100];
    vertices = findVerticesNewMaterialCircle(initPositions.LptDlocs(1,:),settings.polygonSides,1,2);
    initPositions.lpsVertices{1}(:,:) = vertices;
    
    % MAKE IT DISPLAY THE CIRCLE WITH THE SAME AREA AT THE END

    % run the model with given settings

    mainModel(1,settings,initPositions);

end
