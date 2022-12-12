%This is the code to generate market data 
clear;
clc;
myFileName = 'Market.xlsx';
T = 9;
D = [zeros(1,3),300, 300,300,zeros(1,3)];
H = 30 + rand(1,T)*5;
S = 500 + rand(1,T)*100;

xlswrite(myFileName,[D;H;S]);