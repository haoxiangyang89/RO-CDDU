clear;
clc;

ProblemSize = [200,400,600,800,1000,1200];
FileNum = length(ProblemSize);
GridNum = 11;
S = 5000;

PercentA = zeros(FileNum,GridNum);
PercentB = zeros(FileNum,GridNum);
PercentC = zeros(FileNum,GridNum);
for i = 5:FileNum
    myFileName = strcat('E:\F\University of Chicago\Research Chicago\Andy Sun\Matlab Model Code\ErrorTolerance30\myExperiment30\D1N',num2str(ProblemSize(i)),'\');
    for j = 1:GridNum
        [i,j]
        inputFileName = strcat(myFileName,'TAUupper',num2str(j-1),'.mat');
        load(inputFileName);
        
        sumA = sum(sum(Imrp1(1:334,:)));
        sumB = sum(sum(Imrp1(335:668,:)));
        sumC = sum(sum(Imrp1(669:1002,:)));
        
        PercentA(i,j) = sumA/(sumA+sumB+sumC);
        PercentB(i,j) = sumB/(sumA+sumB+sumC);
        PercentC(i,j) = sumC/(sumA+sumB+sumC);
    end
end
plot(0:0.01:0.1,PercentA(5,:),'-*');
hold on;
plot(0:0.01:0.1,PercentB(5,:),'--');
hold on
plot(0:0.01:0.1,PercentC(5,:),'-o');
legend('Type A Resource','Type B Resource','Type C Resource');
xlabel('Budget Control Gamma(t)');
ylabel('Percent of Commitment Amount');
title('Resource Allocation');
