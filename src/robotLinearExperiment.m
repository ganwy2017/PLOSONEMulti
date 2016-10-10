%Experiments with a linear interference and real (with a position) robots
%Delivery point at (0,0)
% Function for experiments with soft deadline
%Nt: number of tasks
%Nr=x*Nt where x=2:2:30
%nExperiments: number of experiment to carry out
function robotLinearExperiment(Nt,nExperiments)
for myR=2:2:30
for t=1:nExperiments
myR
t=38;
fname=sprintf('environments/Environment_SDL_%d',t);
load(fname);


Nr=Nt*myR;

posObj=posObjE(1:Nt,:);


typeObj = typeObjE(1:Nt);

minL = minLE;
maxL = maxLE;
L=LE(1:Nt);

%Utility=unifrnd(100,200,Nt,1);
Utility = UtilityE(1:Nt);

%minCR=1;
%maxCR=14;
%CR=unifrnd(minCR,maxCR,Nr,5); %Robots' load capacity depending on the type
minCR = minCRE;
maxCR = maxCRE;
CR = CRE(1:Nr,1:5);

%vel=repmat(3,Nr,1); %Robots' velocity (same vel. for all robots=3)
%minVel=1;
%maxVel=5;
%vel=unifrnd(minVel,maxVel,Nr,1);
minVel = minVelE;
maxVel = maxVelE;
vel = velE(1:Nr);


%coeficients of the interference function, (distance, a, b) I=ax+b
MaxDistances = 7;
increaseInterference = 1.0;
Coef(1,:) = [137.6065, 10.3428, -15.3982];
Coef(2,:) = [184.0411, 6.1682, -9.2821];
Coef(3,:) = [246.667, 3.4728, -5.3911];
Coef(4,:) = [281.9771, 2.8746, -4.4196];
Coef(5,:) = [328.667, 1.9278, -3.0179];
Coef(6,:) = [359.4595, 1.7190, -2.7518];
Coef(7,:) = [399.1302, 1.3657, -2.0857];
Coef=[Coef(:,1),increaseInterference .* Coef(:,2:3)];

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
        CRij(i,j) = CR(i,typeObj(j)); %For Greedy experiments
         %Get the shortest distance
        [M,I]=min(abs(Coef(:,1)-distance));
        %The interference function is 1000*I
        aij = Coef(I,2)/1000;
        bij = Coef(I,3)/1000;
        F = [F, -1000*Utility(j)*(Kij(i,j) - aij)/L(j)];
    end
    avgDistance = avgDistance / Nt;
end

A=zeros(Nr,Nr*Nt);
AAux=repmat(1,Nt,1)';
 for i=0:Nr-1
     A(i+1,i*Nt+1:(i+1)*Nt) = AAux;
 end
 
b(1:Nr)=1;
b = b';
tic;
X=bintprog(F',A,b); %Remove for debuging 1
myTimes(t,1)=toc;

%Calcular la utilitat total
totalUtilLinear = F*X/(-1000); % Remove for debuging %Alerta falta restar bj/Lj per cada tasca
%totalUtilLinear = 0;
totalU(t,1) = totalUtilLinear;
%------------------------------------------------------------------
%--------------------------Auctions--------------------------------
%------------------------------------------------------------------
avgVel = minVel + ((maxVel-minVel)/2);
avgCR = minCR + ((maxCR - minCR)/2);
avgCR = 0.5 * (avgCR*avgVel/(avgCR*avgVel + avgDistance));
avgL = minL  + ((maxL-minL)/2);
THAuction=avgL/(avgCR * Nt/Nr); %Average time to finish a task
VigParam=[1.0 0.8 0.6 0.4 0.2 0];
nVigParam=length(VigParam);
for pp=0:1
  %VigParam = 0.8;
  %totalU(t,pp+2) = robotAuction(Nt, Nr, posObj, L, Utility, CR, vel, Kij, typeObj, THAuction*0.0, pp, Coef);

  for iVig=1:nVigParam
    tic;
    if (pp ~= 1)
       [totalU(t,(nVigParam*pp)+iVig+1),XNew] = robotAuctionLinear(Nt, Nr, posObj, L, Utility, CR, vel, Kij, typeObj, THAuction*0.0, pp, Coef, VigParam(iVig));
    else
        totalU(t,(nVigParam*pp)+iVig+1) = robotGreedyLoad(Nr,Nt,posObj,L,Utility, Coef,CRij,[]);
    end
    myTimes(t,(nVigParam*pp)+iVig+1)=toc;
  end
  %totalU(t,pp+2) = robotAuctionHDL(Nt, Nr, posObj, L, Utility, CR, vel, Kij, typeObj, THAuction*0.0, pp,[]);
end


%Save partial results if the system can not carry out all 
end
totals(myR/2,:) = mean(totalU);
totalsDev(myR/2,:) = std(totalU);
totalsMax(myR/2,:) = max(totalU);
totalsMin(myR/2,:) = min(totalU);

totalsTimes(myR/2,:) = mean(myTimes);
totalsTimesDev(myR/2,:) = std(myTimes);
totalsTimesMax(myR/2,:) = max(myTimes);
totalsTimesMin(myR/2,:) = min(myTimes);

fname=sprintf('totalU%d_I1_P_%d',Nt, myR);
save(fname,'totalU');
fname=sprintf('totalTimeU%d_I1_P_%d',Nt, myR);
save(fname,'myTimes');

clear totalU;
clear mytimes;
end
