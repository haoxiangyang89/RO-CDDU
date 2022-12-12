# generate DR data based on USC data
using DelimitedFiles, JLD, HDF5;
include("def.jl");

# fileList = ["/Market","/Resource","/Resource2","/Resource3"];
rootPath = "/Users/haoxiangyang/Dropbox/Research_Documents/AP/Collaboration/R-CDDU/julia_code/USC";
resourcePath = rootPath * "/Resource_concave_cost.csv";
data_m = Dict();
data_m["Resource"] = readdlm(resourcePath,',');

T = 9;
n1 = 115;

# need to generate the following field from
cc = data_m["Resource"][1,1:n1];
RD = data_m["Resource"][2,1:n1];
RU = data_m["Resource"][3,1:n1];
pmax = data_m["Resource"][4,1:n1];
pmin = data_m["Resource"][5,1:n1];
upmin = data_m["Resource"][6,1:n1];
downmin = data_m["Resource"][7,1:n1];
alpha = data_m["Resource"][8,1:n1];
beta = data_m["Resource"][9,1:n1];

data = Dict();
for ditem in ["0","300","1400","2000"]
    # market data
    marketPath = rootPath * "/Market_" * ditem * ".csv";
    data_m["Market"] = readdlm(marketPath,',');
    D = data_m["Market"][1,1:T];
    h = data_m["Market"][2,1:T];
    s = data_m["Market"][3,1:T];

    mData = DR_data(T, D, h, s, cc, RD, RU, pmax, pmin, upmin, downmin, alpha, beta);
    data[ditem] = mData;
end
save("testData_USC_concave.jld","mData",data);
