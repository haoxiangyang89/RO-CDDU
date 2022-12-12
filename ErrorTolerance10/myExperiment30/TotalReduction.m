clear;
clc;

ProblemSize = [200,400,600,800,1000,1200];
FileNum = length(ProblemSize);
GridNum = 1;
S = 5000;

PercentA = zeros(FileNum,GridNum);
PercentB = zeros(FileNum,GridNum);
PercentC = zeros(FileNum,GridNum);
for i = 5:FileNum-1
    myFileName = strcat('E:\F\University of Chicago\Research Chicago\Andy Sun\Matlab Model Code\ErrorTolerance10\myExperiment30\D1N',num2str(ProblemSize(i)),'\');
    for j = 1:GridNum
        [i,j]
        inputFileName = strcat(myFileName,'TAUupper',num2str(10),'.mat');
        load(inputFileName);
        
        ScheduleReductionRob = sum(Imrp1)
        
        Demand = myData.('D');
    end
    myFileName = strcat('E:\F\University of Chicago\Research Chicago\Andy Sun\Matlab Model Code\ErrorTolerance10\myExperiment30\D1N',num2str(ProblemSize(i)),'\');
    inputFileName = strcat(myFileName,'TAUupper',num2str(0),'.mat');
    load(inputFileName);
    ScheduleReductionDet = sum(Imrp1);
end
plot(1:length(Demand),Demand,'--');
hold on;
plot(1:length(Demand),ScheduleReductionRob,'-*');
hold on
plot(1:length(Demand),ScheduleReductionDet,'-');
legend('Demand','Robust(0.1)','Deterministic');
xlabel('Time Horizon');
ylabel('Reduction Amount');
title('Total Reduction');
