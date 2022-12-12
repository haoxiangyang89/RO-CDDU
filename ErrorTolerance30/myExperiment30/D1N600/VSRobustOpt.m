clear;
clc;

ProblemSize = [50];
FileNum = length(ProblemSize);
GridNum = 11;
S = 5000;

TimeRecord = zeros(FileNum,GridNum);
SimulatedProfit = zeros(FileNum,GridNum);
for i = 1:FileNum
    myFileName = strcat('E:\F\University of Chicago\Research Chicago\Andy Sun\Matlab Model Code\ErrorTolerance10\myExperiment80\D1N',num2str(ProblemSize(i)),'\');
    for j = 2:GridNum
        [i,j]
        inputFileName = strcat(myFileName,'TAUupper',num2str(j-1),'.mat');
        load(inputFileName);
        
        Zeta = 0.01*(j-1)*sum(myData.('pmax'));
        TimeRecord(i,j) = ImprUppTime;
        N = length(myData.('RD'));
        T = length(myData.('D'));
        [robObjWarm,robP,myRunTime,myGapWarm] = GurobiRobWarm(myData,N,T,Zeta,solvec);
        [ColdrobObj,myRunTimeCold,myGapCold] = GurobiRobCold(myData,N,T,Zeta,myGapWarm);
    
    end
end
