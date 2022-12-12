clear;
clc;
addpath('/apps/gurobi/linux64/matlab');
addpath('E:\Program Files\CPLEX_Optimizer\cplex\matlab\x86_win32\@Cplex')
%Read the data out of excel and populate the arrays in matlab
M=3;
T = 9; 
N = 334;
myData = struct();


%Read mmarket data
myDataPath = 'Market';
myData1 = xlsread(myDataPath);

myData.('D') = myData1(1,1:T);
myData.('HoldingCost') = myData1(2,1:T);
myData.('ShortageCost') = myData1(3,1:T);

%Read resource data
myDataPath = strcat('Resource',num2str(1),'.xlsx');
myData1 = xlsread(myDataPath);
myData.('CommitmentCost') = myData1(1,1:N);
myData.('RD') = myData1(2,1:N);
myData.('RU') = myData1(3,1:N);
myData.('pmax') = myData1(4,1:N);
myData.('pmin') = myData1(5,1:N);
myData.('Upupmin') = myData1(6,1:N);
myData.('Downdownmin') = myData1(7,1:N);
myData.('alpha') = myData1(8,1:N);
myData.('beta') = myData1(9,1:N);
for i = 2:M 
    myDataPath = strcat('Resource',num2str(i),'.xlsx');
    myData1 = xlsread(myDataPath);
    myData.('CommitmentCost') = [myData.('CommitmentCost'),myData1(1,1:N)];
    myData.('RD') = [myData.('RD'),myData1(2,1:N)];
    myData.('RU') = [myData.('RU'),myData1(3,1:N)];
    myData.('pmax') = [myData.('pmax'),myData1(4,1:N)];
    myData.('pmin') = [myData.('pmin'),myData1(5,1:N)];
    myData.('Upupmin') = [myData.('Upupmin'),myData1(6,1:N)];
    myData.('Downdownmin') = [myData.('Downdownmin'),myData1(7,1:N)];
    myData.('alpha') = [myData.('alpha'),myData1(8,1:N)];
    myData.('beta') = [myData.('beta'),myData1(9,1:N)];
end
%Discussion 1
%Deterministic solution

[detObj, detP,detTime,solvec] = GurobiDet(myData,3*N,T);

% tic
% [detObj, detP,detTime] = GurobiDetBackup(myData,3*N,T);
% toc
% myTitle = strcat('TAUupper',num2str(0));
% save(myTitle,'detP','myData','detTime');

for i = 0:10
%Set the Zeta value 
Zeta = 0.01*i*sum(myData.('pmax'));
%Discussion 2
%Improved upperbound solution for robust optimization model initial on p
Imrp1=detP;
error = 1;
provObj = 0;
curImprObj = 10;
ImprUppTime =0;
prevSol = solvec;
while error>10^(-3)
    [obj,Pi,Lam,Mu,Tht,myRunTime]=ImprDual(myData,3*N,T,Imrp1,Zeta);
    ImprUppTime = ImprUppTime+myRunTime;
    prevObj = curImprObj;
    [curImprObj,Imrp1,solvec,myRunTime]=ImprMst(myData,3*N,T,Pi,Lam,Mu,Tht,Zeta,prevSol);
    prevSol = solvec;
    ImprUppTime = ImprUppTime+myRunTime;
    error = abs(curImprObj - prevObj)/abs(prevObj);
    %[Y,hahaY1,hahaY2] = RobFeasibilityChecking(myData,3*N,T,Pi,Lam,Mu,Tht,Zeta,p)
end

%Discussion 3
%Improved upperbound solution for robust optimization model initial on pi 

% error = 1;
% provObj = 0;
% curImprObj2 = 10;
% for t = 1:T
%     Mu(t) = 0;
%     Tht(t) = 0; 
%     for i = 1:3*N
%         Pi(i,t) = abs(myData.('HoldingCost')(t) - myData.('CommitmentCost')(i));
%         Lam(i,t) = abs(myData.('ShortageCost')(t) + myData.('CommitmentCost')(i));
%     end
% end

% tic
% while error>10^(-3)
%     prevObj = curImprObj2;
%     [curImprObj2,Imrp2,solvec]=ImprMst(myData,3*N,T,Pi,Lam,Mu,Tht,Zeta);
%     
%     [obj,Pi,Lam,Mu,Tht]=ImprDual(myData,3*N,T,Imrp2,Zeta);
%     error = abs(curImprObj2 - prevObj)/abs(prevObj);
%     %[Y,hahaY1,hahaY2] = RobFeasibilityChecking(myData,3*N,T,Pi,Lam,Mu,Tht,Zeta,p)
% end
% ImprUppTime2 = toc

%Discussion 4
%Upperbound solution for Robust Optimization Model
% Upp=detP;
% 
% error=1;
% prevObj = 0; 
% curObj = 10;
% tic
% while error > 10^(-3)
%     [obj,U,V,R,S,W,Up,Vp,Rp,Sp,Wp]=SubProblemABRS(myData,3*N,T,Upp,Zeta);
%     prevObj = curObj;
%     [curObj,Upp]=SubProblemP(myData,3*N,T,U,V,R,S,W,Up,Vp,Rp,Sp,Wp,Zeta);
%     error = abs(curObj - prevObj)/abs(prevObj);
% end
% NormUppTime = toc


%Discussion 5
%Global Solution for Robust Optimization Model
%RobFeasibilityChecking(myData,3*N,T,Zeta,x);
%We feed the upperbound solution to the global model 
myTitle = strcat('TAUupper',num2str(i));
save(myTitle,'Imrp1','myData','ImprUppTime','solvec');
% [robObj,robP,myRunTime] = GurobiRobWarm(myData,3*N,T,Zeta,solvec);
% robTimeWarmStart = myRunTime+ImprUppTime;
% myTitle = strcat('TAURobWarm',num2str(i));
% save(myTitle,'robP','myData','robTimeWarmStart');
% 
% [robObj,robP,myRunTime] = GurobiRobCold(myData,3*N,T,Zeta);
% robTimeColdStart = myRunTime;
end

