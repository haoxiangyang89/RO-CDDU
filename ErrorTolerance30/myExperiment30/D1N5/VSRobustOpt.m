clear;
clc;

ProblemSize = [50];
FileNum = length(ProblemSize);
GridNum = 11;
S = 5000;

BaseTime = [3.04]
TimeRecord = zeros(FileNum,GridNum);
TimeRecordRob = zeros(FileNum,GridNum);
for i = 1:FileNum
    myFileName = strcat('E:\F\University of Chicago\Research Chicago\Andy Sun\Matlab Model Code\ErrorTolerance30\myExperiment30\D1N',num2str(ProblemSize(i)),'\');
    for j = 11:GridNum
        [i,j]
        inputFileName = strcat(myFileName,'TAUupper',num2str(j-1),'.mat');
        load(inputFileName);
        TimeRecord(i,j) = BaseTime(i)+ImprUppTime;  
        
        Zeta = 0.01*(j-1)*sum(myData.('pmax'));
        TimeRecord(i,j) = ImprUppTime;
        N = length(myData.('RD'));
        T = length(myData.('D'));
        
        robObjWarm=0;
        %[robObjWarm,robP,myRunTime,myGapWarm] = GurobiRobWarm(myData,N,T,Zeta,solvec);
        [ColdrobObj,myRunTimeCold,myGapCold] = GurobiRobCold(myData,N,T,Zeta,robObjWarm);
        TimeRecord2(i,j) = myRunTimeCold;
    end
end
