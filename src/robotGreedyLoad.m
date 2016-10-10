%[totalU]=robotGreedy(Nr,Nt,posObj,L,Utility,Kij,Coef,DL)
%Greedy Experiments where each robot select the best task for it
%Nr: number of robots
%Nt: number of tasks
%posObj= Ntx2 matrix with the positions of the objects
%L: vector of Nt element with the load of each object
%Utility: vector of Nt elements with the utility of each task
%Coef: interference coefficients. it follows the same format as in robotLinearExperiment function
%CCij: NrxNt matrix with the capacity of each robot to execute each task
%DL: Nt vector, each component is the deadline of a task. If the vector is empty the tasks are soft-deadline.

%totalU: total obtained utility 
function [totalU]=robotGreedyLoad(Nr,Nt,posObj,L,Utility,Coef,CCij,DL)
   [maxU, maxJ] = max(CCij');
   totalU = 0;
   totalRobotsG = 0;
   AssignedRobots=[];
   for j=1:Nt
       robotsTask=find(maxJ==j); %index dels robots que feran la tasca j
       AssignedRobots = [AssignedRobots, robotsTask];
       utilitiesTask = CCij(robotsTask,j);
       Cj = sum(utilitiesTask);
       
       
       %Get total Utility
       %1. Get the interference for num. robots = nRobots
       n = length(robotsTask);
       
       if (n == 0) 
           continue;
       end
        distance = sqrt((posObj(j,1).^2)+(posObj(j,2).^2));
        [M,I]=min(abs(Coef(:,1)-distance));
        %La funcio de interferencia esta multiplicada per 1000
        aij = Coef(I,2)/1000;
        bij = Coef(I,3)/1000;
        %interference=polyval([aij bij], n);
        %if (n==1) 
        %    interference = 0;
        %end
        

        totalRobotsG = totalRobotsG + n;
        interference=polyval([aij,0], n);
        Cj = Cj - interference;
        predictedT = L(j)/Cj;
        
        if (predictedT < 0)
           continue;
        end
        
        if isempty(DL) %SDL experiment
            totalU = totalU + Utility(j)/predictedT;
        else %HDL
             if (predictedT <= DL(j))
                totalU = totalU + Utility(j);
                %disp('okiDL');
                %pause;
             end
        end
   end
   AssignedRobotsAux = unique(AssignedRobots);
   if length(AssignedRobotsAux) < length(AssignedRobots)
      pause;
   end
   
end
