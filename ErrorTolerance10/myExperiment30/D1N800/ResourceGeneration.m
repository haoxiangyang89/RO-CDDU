clear;
clc;
%This is the code to generate resource data 
%Number of resources in each type 
N = 267;


OUnc = [0.5,0.3,0.1];
%Data source 1
A = zeros(9,N);
%Commitment Cost
A(1,:) = round(100*rand(1,N)*2)/100+22;
%RD
A(2,:) = round(rand(1,N)*2)+5;
%RU
A(3,:) = round(rand(1,N)*2)+5;
%pmax
A(4,:) = round(rand(1,N)*5)+15;
%pmin
A(5,:) = round(rand(1,N))+4;
%Upupmin
A(6,:) = round(rand(1,N)*2)+2;
%Downdownmin
A(7,:) = round(rand(1,N)*2)+2;
%alpha
A(8,:) = ones(1,N)*(-OUnc(1));
%beta
A(9,:) = ones(1,N)*(OUnc(1));
myFileName = strcat('Resource',num2str(1),'.xlsx');
xlswrite(myFileName, A);

%Data source 2
A = zeros(9,N);
%Commitment Cost
A(1,:) = round(100*rand(1,N)*2)/100+20;
%RD
A(2,:) = round(rand(1,N)*2)+5;
%RU
A(3,:) = round(rand(1,N)*2)+5;
%pmax
A(4,:) = round(rand(1,N)*5)+15;
%pmin
A(5,:) = round(rand(1,N))+4;
%Upupmin
A(6,:) = round(rand(1,N)*2)+2;
%Downdownmin
A(7,:) = round(rand(1,N)*2)+2;
%alpha
A(8,:) = ones(1,N)*(-OUnc(2));
%beta
A(9,:) = ones(1,N)*(OUnc(2));

myFileName = strcat('Resource',num2str(2),'.xlsx');
xlswrite(myFileName, A);

%Data source 3
A = zeros(9,N);
%Commitment Cost
A(1,:) = round(100*rand(1,N)*2)/100+18;
%RD
A(2,:) = round(rand(1,N)*2)+5;
%RU
A(3,:) = round(rand(1,N)*2)+5;
%pmax
A(4,:) = round(rand(1,N)*5)+15;
%pmin
A(5,:) = round(rand(1,N))+4;
%Upupmin
A(6,:) = round(rand(1,N)*2)+2;
%Downdownmin
A(7,:) = round(rand(1,N)*2)+2;
%alpha
A(8,:) = ones(1,N)*(-OUnc(3));
%beta
A(9,:) = ones(1,N)*(OUnc(3));

myFileName = strcat('Resource',num2str(3),'.xlsx');
xlswrite(myFileName, A);


