%This is the code to generate market data 
clear;
clc;
myFileName = 'Market.xlsx';
T = 9;
D = [zeros(1,3),4800, 4800,4800,zeros(1,3)];
H = 30 + rand(1,T)*5;
S = 1000 + rand(1,T)*500;

xlswrite(myFileName,[D;H;S]);