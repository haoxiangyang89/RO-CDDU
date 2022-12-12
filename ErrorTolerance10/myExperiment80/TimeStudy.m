clear;
clc;

ProblemSize = [200,400,600,800,1000,1200];
FileNum = length(ProblemSize);
GridNum = 11;
S = 5000;
BaseTime = [47.40,120.52,205.26,296.16,974.61,952.02,1002.36]
ComputationTime = zeros(FileNum,GridNum);
for i = 1:FileNum
    myFileName = strcat('E:\F\University of Chicago\Research Chicago\Andy Sun\Matlab Model Code\ErrorTolerance10\myExperiment80\D1N',num2str(ProblemSize(i)),'\');
    for j = 1:GridNum
        [i,j]
        inputFileName = strcat(myFileName,'TAUupper',num2str(j-1),'.mat');
        load(inputFileName);
        ComputationTime(i,j) = BaseTime(i)+ImprUppTime;    
    end
end
mesh(ComputationTime)
