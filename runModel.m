function runModel(m)
% m is number of iterations

if nargin < 1
    m = 50;
end

% set the random seed
%rng(27)

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
settings.splitFlag = 0;
settings.plotEvery = 0;

settings.insRateProtein = 1500;
settings.insRateLPS = 1500;
settings.insRateBAM = 10;
settings.insRateLptD = 0.022;  % = 300/15(30^3) p46 lab book 4

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

%initPositions.BAMlocs = [[100,100]];

initPositions.BAMlocs = [];
initPositions.proteinVertices = {};
initPositions.LptDlocs = [];
initPositions.lpsVertices = {};

for i = 1500
    i
    % set up the settings that vary
    settings.insRateProtein = i;
    
    settings.saveLocation = [bigSaveFolder,'/test-proteinInsertionRate',num2str(i)];

    % run the model with given settings
    
    %try
    model = mainModel(0,settings,initPositions);
    visualiseSimple(model.BAMlocs,model.proteinVertices,model.LptDlocs,model.lpsVertices,model.rightEdge,model.leftEdge);
    %catch e
%         e
%         visualiseSimple(model.BAMlocs,model.proteinVertices,model.LptDlocs,model.lpsVertices,model.rightEdge,model.leftEdge);
%     end
end
