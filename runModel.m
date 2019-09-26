function model = runModel(m,useGivenInsertions)
% m is number of iterations
% useGivenInsertions is 0 if not and a file location if yes

if nargin < 1
    m = 3000;
end

if useGivenInsertions ~= 0
    data = load(useGivenInsertions);
    oldModel = data.model;
    
    initPositions.BAMlocs = oldModel.BAMlocs;
    
    initPositions.LptDlocs = oldModel.LptDlocs;
else
    initPositions.BAMlocs = [];
    
    initPositions.LptDlocs = [];
end

% set the random seed and record it
rng('shuffle'); % first re-set the random generator, otherwise get same results each time
randSeed = randi([0,1000],1); % pick a random seed
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
settings.plotEvery = 200;

  % 0.022 = 300/15(30^3) p46 lab book 4

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

settings.initialMembraneArea = settings.membraneCircumference*settings.currentMaxLen;
settings.doublingTime = 30;
settings.areaFractionProtein = 0.3;
settings.areaFractionLPS = 1-settings.areaFractionProtein;
settings.BAMAverage = 445;

settings.LptDaverage = 100:50:300;

for i = 1:numel(settings.LptDaverage)
    i
    % calculate settings
    L = settings.LptDaverage(i)

    QL = settings.areaFractionLPS*settings.initialMembraneArea/(settings.doublingTime*L);
    kL = 9*L/(13*settings.BAMAverage*settings.doublingTime);

    NB0 = 6*settings.areaFractionLPS*settings.initialMembraneArea/(13*(settings.doublingTime^2)*kL*QL);
    kB = 6*settings.areaFractionLPS*settings.initialMembraneArea/(13*(settings.doublingTime^3)*kL*QL);

    NL0 = 9*settings.areaFractionLPS*settings.initialMembraneArea/(13*settings.doublingTime*QL);

    Qp = (13*settings.areaFractionProtein*settings.doublingTime*kL*QL)./(9*settings.areaFractionLPS);
    
    settings.insRateProtein = Qp;
    settings.insRateLPS = QL;
    settings.insRateBAM = kB;
    settings.insRateLptD = kL;
    
    [BAMlocs,LptDlocs] = makeInitialBAMlocsLptDlocs(settings.membraneCircumference,NB0,NL0);
    
    initPositions.BAMlocs = BAMlocs;
    initPositions.LptDlocs = LptDlocs;
    
    settings.saveLocation = [bigSaveFolder,'/test-klptD',num2str(kL),'-Qlps-',num2str(QL)];

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

end

function [BAMlocs,LptDlocs] = makeInitialBAMlocsLptDlocs(membraneCircumference,NBAM,NLptD)
    BAMlocs = zeros(NBAM,2);
    for i=1:NBAM
        newBAMloc = rand(1,2);
        newBAMloc(1) = 2000*newBAMloc(1)-1000;
        newBAMloc(2) = membraneCircumference*newBAMloc(2);
        BAMlocs(i,:) = newBAMloc;
    end
    
    LptDlocs = zeros(NLptD,2);
    membraneCircumference = 500*pi;
    for i=1:NLptD
        newLptDloc = rand(1,2);
        newLptDloc(1) = 2000*newLptDloc(1)-1000;
        newLptDloc(2) = membraneCircumference*newLptDloc(2);
        LptDlocs(i,:) = newLptDloc;
    end
end
