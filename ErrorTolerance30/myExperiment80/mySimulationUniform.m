clear;
clc;

ProblemSize = [200,400,600,800,1000,1200];
FileNum = length(ProblemSize);
GridNum = 11;
S = 5000;


SimulatedProfit = zeros(FileNum,GridNum);
for i = 5:FileNum
    myFileName = strcat('E:\F\University of Chicago\Research Chicago\Andy Sun\Matlab Model Code\ErrorTolerance30\myExperiment80\D1N',num2str(ProblemSize(i)),'\');
    for j = 1:GridNum
        [i,j]
        inputFileName = strcat(myFileName,'TAUupper',num2str(j-1),'.mat');
        load(inputFileName);
        N = length(myData.('RD'));
        T = length(myData.('D'));
        
        myTempProfit = 0;
        for k = 1:S
            LAlpha = myData.('alpha')'*ones(1,T);
            UBeta = myData.('beta')'*ones(1,T);
            RealizedReduction = sum((1+(rand(N,1)*ones(1,T).*(UBeta-LAlpha)+LAlpha)).*Imrp1,1); 
            Short = max(myData.('D') - RealizedReduction,0);
            More = max(RealizedReduction- myData.('D'),0); 
            ActualProfit = -myData.('ShortageCost')*Short'-myData.('HoldingCost')*More'...
                +myData.('CommitmentCost')*sum(Imrp1,2);
            
            myTempProfit = myTempProfit+ActualProfit;
        end
        SimulatedProfit(i,j) = myTempProfit/S;
    end
end
plot(0:0.01:0.1,SimulatedProfit(5,:))
xlabel('Budget Control Gamma(t)');
ylabel('Simulated Profit');
title('Simulated Profit vs Budget Control');