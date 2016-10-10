
%[auctionTotalU]=robotAuctionHDL(Nt, Nr, posObj, L, Utility, CR, vel, Kij, typeObj, THAuction, auctionType, DL, Coef, VigParam)
%Auction experiments with iterative double round auction process and soft-deadlines
%Nt: numbero of tasks
%Nr: number of robots
%posObj: objects' positions (Nt x 2)
%L: Objects' weights (Nt)
%Utility: tasks' utilities (Nt
%CR: robots' load capacity  (Nr)
%vel: velocity of the robots (Nr)
%Kij: K(ij)= ith robot's work capacity for task jth. (NrxNt)
%typeObj: Nt vector. typeObj(i): type of the ith Object.
%THAuction: overhead over the last best utility. Set 0 for PlosOne
%experiments
%auctionType: 0 if the robots follows MDRA o SDRA algorithm, otherwise random selection.
%Coef: coeficients of interference function I.
%VigParam: lambda_B parameter of SRDA and MDRA algortihms

%auctionTotalU: total obtained utility


function [auctionTotalU]=robotAuctionHDL(Nt, Nr, posObj, L, Utility, CR, vel, Kij, typeObj, THAuction, auctionType, DL, Coef, VigParam)

    function [predictedU predictedT totalCapacityNew DLOK] = predictedUtility(totalCapacity, n, robotIndex, taskIndex)
        %coeficients of the interference function, (distance, a, b) I=ax+b
        DLOK = 0;
        MaxDistances = 7;
        %Coef(1,:) = [137.6065, 10.3428, -15.3982];
        %Coef(2,:) = [184.0411, 6.1682, -9.2821];
        %Coef(3,:) = [246.667, 3.4728, -5.3911];
        %Coef(4,:) = [281.9771, 2.8746, -4.4196];
        %Coef(5,:) = [328.667, 1.9278, -3.0179];
        %Coef(6,:) = [359.4595, 1.7190, -2.7518];
        %Coef(7,:) = [399.1302, 1.3657, -2.0857];
        
        distance = sqrt((posObj(j,1).^2)+(posObj(j,2).^2));
        [M,I]=min(abs(Coef(:,1)-distance));
        %La funcio de interferencia esta multiplicada per 1000
        aij = Coef(I,2)/1000;
        bij = Coef(I,3)/1000;
        %interference=polyval([aij bij], n);
        %if (n==1) 
        %    interference = 0;
        %end
        interference=polyval([aij], n);
        totalCapacity;
        totalCapacityNew = totalCapacity + (0.5 * (CR(robotIndex,typeObj(j))*vel(robotIndex)/((vel(robotIndex)*CR(i,typeObj(j))) + distance)));
        predictedT = L(taskIndex)/(totalCapacityNew-interference);
        if (length(DL) > 0)
            %fprintf('checking DL: %g, index=%i, predicted=%g\n',DL(taskIndex), taskIndex, predictedT);
            if (predictedT <= DL(taskIndex))
                DLOK = 1;
                predictedU = Utility(taskIndex);
            else
                predictedU = 0;
            end
        else
           predictedU = Utility(taskIndex)/predictedT;  
        end
    end

 function [totalUtility] = getTotalUtility(taskPerRobot, leadersPerTask)
        totalUtility = 0;
        for j=1:Nt
            totalUtilityT=0;
            %Get the leader
            leaderIndex = leadersPerTask(j);
            nRobots = 1;
            [predictedU predictedT totalCapacity DLOK] =  predictedUtility(0,nRobots,leaderIndex,j);
            for i=1:Nr
                if (taskPerRobot(i) ==j)
                    
                    nRobots = nRobots + 1;
                    [predictedU predictedT totalCapacity DLOK] =  predictedUtility(totalCapacity,nRobots,i,j);
                end
            end
            totalUtilityT = totalUtilityT + predictedU;
            
            if (length(DL) > 0) %Hard-deadline function
            totalUtilityT = 0;
              if (DL(j)>=predictedT)
                totalUtilityT = Utility(j);
              end
            end
            totalUtility = totalUtility + totalUtilityT;
        end
    end

   %fprintf('Starting Auction...\n');
   overheadTH = 1.1;
   %1st. Select as a leader the robot. The robot with the max capacity
   auxIndex = [1:Nr]';
   KijAux = [Kij auxIndex];
   leaders = [];
   for j=1:Nt
      [maxRob leaderI] = max(KijAux(:,j));
      leaders = [leaders, KijAux(leaderI,Nt+1)];
      KijAux(leaderI,:)=[];
   end
   
   
   %2nd. each leader tries to create each own working group
   %All robots send its capacity for that task and the leaders select the
   %best team until a TH is reached
   
   %bids(j,i): robot's index for task j
   bids = repmat(-1,Nt,Nr);
   for j=1:Nt
       insideIndex=1;
       for i=1:Nr
           if (length(find(leaders==i))==0)
             %El robot no es leader
             bids(j,insideIndex) = i;
             insideIndex = insideIndex + 1;
           end
       end
   end
   
   %3th. Selecting the best robots
   %myleadrs(Nr) 
   %myLeaders(i).list. list of canditates to become ith's robot leader
   %myLeaders=[Nr,Nt]; for each robot it contains the index of the task!!!!
   myLeaders = repmat(-1,Nr,Nt);
   myLeadersUtil = repmat(-1,Nr,Nt);
   nLeadersRobot = repmat(0,1,Nr);
   for j=1:Nt
      
       listIndex = bids(j,:);
       i=1;
       while ((i<=Nr) && (listIndex(i)>=1))
           listKij(i) = Kij(listIndex(i),j);
           i = i +1;
       end
       [value, index] = sort(listKij,2,'descend');
       
       %Get the total capacity
       leaderIndex = leaders(j);
       totalCapacity = 0;
       [newUtility newTime totalCapacity, DLOK] = predictedUtility(totalCapacity,1,leaderIndex,j);
       oldUtility = newUtility;
       n=1;
       %if (DLOK==1)
       %    fprintf('Groups limit reached 1\n');
       %end
       while ((n<=length(index)) && (DLOK == 0))
          %fprintf('while 1');
          if (auctionType == 2)
              [newUtility newTime totalCapacity, DLOK] = predictedUtility(totalCapacity,1,index(n),j);
          else
              [newUtility newTime totalCapacity, DLOK] = predictedUtility(totalCapacity,n+1,index(n),j);
          end
          %if (DLOK==1)
              %fprintf('Groups limit reached\n');
          %    break;
          %end
          %List with the members working group of task j
          groups(j,n).idRobot = listIndex(index(n));
          %if (auctionType==2)
          %    groups(j,n).utility = rand();
          %else
             groups(j,n).utility = (newUtility-oldUtility);
          %end
          
          %end
          oldUtility = newUtility;
          n = n+1;
       end
       groups(j,n).idRobot = (-1);
       groups(j,n).utility = (-1);
       
        %fill the list of leaders for each robot
        for k=1:n-1
            if ((auctionType==0) || (auctionType==2))
               groups(j,k).utility = oldUtility;
            end
            nLeadersRobot(groups(j,k).idRobot) = nLeadersRobot(groups(j,k).idRobot) + 1;
            myLeaders(groups(j,k).idRobot, nLeadersRobot(groups(j,k).idRobot)) = j;
            myLeadersUtil(groups(j,k).idRobot, nLeadersRobot(groups(j,k).idRobot)) = groups(j,k).utility;
        end
   end
   %leaders
   %myLeaders
   %groups(j,:) members working group of task j
   %leaders(j) leader of task j
   %Second phase: the robots select the best leaders
   for i=1:Nr
       k=1;
       maxUtilityR = -10000;
       indexTaskMax = (-1);
       %i
       maxUtility=max(myLeadersUtil(i,:));
       while ((k<=Nt) && (myLeaders(i,k)>=1))
           
           if ((myLeadersUtil(i,k)>=(VigParam*maxUtility)) && (Kij(i,myLeaders(i,k)) > maxUtilityR))
               indexTaskMax=myLeaders(i,k);
               maxUtilityR = Kij(i,indexTaskMax);
           end
           k=k+1;
       end
       myLeaders(i,1) = indexTaskMax;
   end
   %myLeadersUtil
   jj=myLeaders(:,1)';
   
   auctionTotalU = getTotalUtility(myLeaders(:,1)', leaders);
   
end  
