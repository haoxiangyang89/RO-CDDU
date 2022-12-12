clear;
clc;

ProblemSize = [200,400,600,800,1000,1200];
FileNum = length(ProblemSize);
GridNum = 11;
S = 5000;
BaseTime = [40.7,154.32,167.22,653.81 ,727.2,1001.90];
TimeRecord = zeros(FileNum,GridNum);
ObjectiveRecrod = zeros(FileNum,GridNum);
for i = 1:FileNum
    myFileName = strcat('E:\F\University of Chicago\Research Chicago\Andy Sun\Matlab Model Code\ErrorTolerance30\myExperiment30\D1N',num2str(ProblemSize(i)),'\');
    for j = [2,5,8]
        [i,j]
        inputFileName = strcat(myFileName,'TAUupper',num2str(j-1),'.mat');
        load(inputFileName);
        
%        Zeta = 0.01*(j-1)*sum(myData.('pmax'));
        TimeRecord(i,j) = ImprUppTime+BaseTime(i);
%         N = length(myData.('RD'));
%         T = length(myData.('D'));
%         [robObjWarm,robP,myRunTime,myGapWarm] = GurobiRobWarm(myData,N,T,Zeta,solvec);
%         [ColdrobObj,myRunTimeCold,myGapCold] = GurobiRobCold(myData,N,T,Zeta,myGapWarm);
%     
    end
end
