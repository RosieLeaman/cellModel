function mainModel

% model settings

insRateProtein = 0;
insRateLPS = 0;

time = 0;
dt = 0.01;
maxTime = 0.05;

% set up the model

BAMlocs = [];
LptDlocs = [];
proteinVertices = [];
lpsVertices = [];

% make an insertion 

vertices = findVerticesNewMaterial([1,1],20,0.1);

% while time is less than max time
figure;hold on;
plot(vertices(:,1),vertices(:,2),'o')
while time < maxTime

    % see if we need to insert a new BAM
    
        % if yes add a BAM to a new randomly chosen location
        
        % add new material (protein)
       
    % see if we need to add new lptD (or how many?)
    
        % if yes pick a BAM
        
        % add new material (lipid)
              
    % move polygon vertices
    %for j=13:13
    for j=1:size(vertices,1)
        %flow = calcFlow(vertices(j,:),[[1,1]],1,[],0,1);
        flow = calcFlow(vertices(j,:),[[0.8,0.8];[1,1]],1,[],0,1);
        vertices(j,:) = findNewVertexPosition(vertices(j,:),flow,dt);
    end
    
    plot(vertices(:,1),vertices(:,2),'o')
    
    
    % move insertion points
     
    
    % resolve boundary issues
    
    
    % calculate anything you want to calculate
    
    
    % increase time
    
    time = time + dt;
    
end
plot([0.8,1],[0.8,1],'x','linewidth',2)
%plot([1],[1],'x','linewidth',2)
legend('1','2','3','4','5','6')
    