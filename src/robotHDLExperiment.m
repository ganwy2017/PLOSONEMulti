%robotHDLExperiment(Nt, Nr, numSimHDL, incInterference)
%Experiments with hard-deadline functions
%Delivery point at (0,0)
%Nt: number of tasks
%Nr: number of robots
%numSimHDL: number of simulations to perform
%incInterference: increment of the orginal interference, equals to 1 for all experiments

function robotHDLExperiment(Nt, Nr, numSimHDL, incInterference)
%for myR=2:2:6
nErrors = 0;
nExcepcions = 0;
for t=1:numSimHDL

    
fname=sprintf('environments/Environment_HDL_%d',t);
load(fname);
Nt %=4%30; %25
Nr %=20 %15%Nt*myR;
nErrors
nExcepcions
t
%posObj=unifrnd(-283,283,Nt,2); %Objects' positions posRobot=100x2
posObj=posObjE(1:Nt,:);

%typeObj=round(4*rand(1,Nt)+1);
typeObj = typeObjE(1:Nt);

%minL = 10;
%maxL = 30;
%L=unifrnd(minL,maxL,Nt,1); %Object's weigth
minL = minLE;
maxL = maxLE;
L=LE(1:Nt);

%Utility=unifrnd(100,200,Nt,1); %repmat(100,1,Nt); %unifrnd(100,200,Nt,1);
Utility = UtilityE(1:Nt);

%minCR=1;
%maxCR=5; %14;
%CR=unifrnd(minCR,maxCR,Nr,5); %Robots' load capacity depending on the type
minCR = minCRE;
maxCR = maxCRE;
CR = CRE(1:Nr,1:5);

% minVel=1;
% maxVel=5;
% vel=unifrnd(minVel,maxVel,Nr,1);
minVel = minVelE;
maxVel = maxVelE;
vel = velE(1:Nr);


%Deadlines
% minDL=50; %34 ;
% maxDL=400;%267;
% DL=unifrnd(minDL,maxDL,Nt,1);

minDL=minDLE;
maxDL=maxDLE;
DL=DLE(1:Nt);

%coefficients of the interference function, (distance, a, b) I=ax+b
MaxDistances = 7;
Coef(1,:) = [137.6065, 10.3428, -15.3982];
Coef(2,:) = [184.0411, 6.1682, -9.2821];
Coef(3,:) = [246.667, 3.4728, -5.3911];
Coef(4,:) = [281.9771, 2.8746, -4.4196];
Coef(5,:) = [328.667, 1.9278, -3.0179];
Coef(6,:) = [359.4595, 1.7190, -2.7518];
Coef(7,:) = [399.1302, 1.3657, -2.0857];
Coef=[Coef(:,1),incInterference .* Coef(:,2:3)];

%min F*Xt
%s.t.
%AX<=b
F=[];
for i=1:Nr
    avgDistance = 0;
    for j=1:Nt
        distance = sqrt((posObj(j,1).^2)+(posObj(j,2).^2));
        avgDistance = avgDistance + distance;
        Kij(i,j) = 0.5 * (CR(i,typeObj(j))*vel(i)/((vel(i)*CR(i,typeObj(j))) + distance));
         %Get the shortest distance
        [M,I]=min(abs(Coef(:,1)-distance));
        %La funcio de interferencia esta multiplicada per 1000
        aij = Coef(I,2)/1000;
        bj(j) = Coef(I,3)/1000;
        Cij(i,j)= Kij(i,j) - aij;
        
    end
    avgDistance = avgDistance / Nt;
end

%Create the optimizator proper structures F, A, b
minM=min((-L./DL)'-bj)
MT = diag(repmat(minM,1,Nt));
CT=[];
for i=1:Nr
    CA=-Cij(i,:); %+aj;
    CT=[CT diag(CA)];
end
CT = [MT CT];
A=zeros(Nr,Nr*Nt);
AAux=repmat(1,Nt,1)';
 for i=0:Nr-1
     A(i+1,i*Nt+1:(i+1)*Nt) = AAux;
 end
 
 A = [zeros(Nr,Nt) A];
 A = [CT;A];
 b=[(-L./DL)'-bj,repmat(1,1,Nr)];
 %F=[repmat(1,1,Nt)./Utility repmat(0,1,Nr*Nt)];
  F=[Utility' repmat(0,1,Nr*Nt)];
 
 newOptions = optimset('MaxTime', 99999, 'MaxIter',900000*(Nt+Nr*Nt), 'MaxRLPIter', 999999*(Nt+Nr*Nt));
 
 X0=repmat(0,1,Nt+Nt*Nr);
 try
   tic;
   [X sortida1 sortida2]=bintprog(F,A,b); %,[],[],X0,newOptions);
 catch
   sortida2 = 2;
   nExcepcions = nExcepcions + 1; 
 end

%Calcular la utilitat total
if (sortida2 == 1)
  myTimes(t,1)=toc;
  totalUtilLinear = sum((1-X(1:Nt))' .* Utility')
  totalU(t,1) = totalUtilLinear;
else
  totalU(t,1) = (-1234.0);
  myTimes = (-1234.0);
  nErrors=nErrors + 1;
  errorsSim(nErrors) = t;
end
%------------------------------------------------------------------
%--------------------------Auctions--------------------------------
%------------------------------------------------------------------

fprintf('Auctions.....\n');
VigParam=[1.0 0.8 0.6 0.4 0.2 0];
nVigParam=length(VigParam);
for pp=0:1
  avgVel = minVel + ((maxVel-minVel)/2);
  avgCR = minCR + ((maxCR - minCR)/2);
  avgCR = 0.5 * (avgCR*avgVel/(avgCR*avgVel + avgDistance));
  avgL=minL + ((maxL-minL)/2);
  THAuction=avgL/(avgCR * Nt/Nr); %Average time to finish a task
  for iVig=1:nVigParam
    tic;
    if (pp ~= 1)
       totalU(t,(nVigParam*pp)+iVig+1) = robotAuctionHDL(Nt, Nr, posObj, L, Utility, CR, vel, Kij, typeObj, THAuction*0.0, pp,DL, Coef, VigParam(iVig));
    else
        totalU(t,(nVigParam*pp)+iVig+1) = robotGreedyLoad(Nr,Nt,posObj,L,Utility,Coef,Kij,DL);
    end
    myTimes(t,(nVigParam*pp)+iVig+1)=toc;
  end
end

clear X0
clear b
clear F
clear X
clear sortida1
%if (sortida2 == 1) 
%  totalU(t,pp+3)=robotVigHDL(Nt, Nr, posObj, L, Utility, CR, vel', Kij, typeObj, Coef, MaxDistances, 1, DL);
%else
%  totalU(t, pp+3) = (-1234.0);
%end

fname=sprintf('totalResHDE_%d_%d_I%d',Nt,Nr,incInterference)
save(fname);
end
end
