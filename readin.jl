# read in data
using DelimitedFiles, JLD, HDF5;
include("def.jl");

firstsub = ["/ErrorTolerance10","/ErrorTolerance30"];
secondsub = ["/myExperiment30","/myExperiment80"];
lastsub = ["/D1N5","/D1N20","/D1N50","/D1N100","/D1N200","/D1N400","/D1N600","/D1N800","/D1N1000","/D1N1200"];
lastn = [5,20,50,100,200,400,600,800,1000,1200];
fileList = ["/Market","/Resource1","/Resource2","/Resource3"];
rootPath = "/Users/haoxiangyang/Dropbox/Research_Documents/AP/Collaboration/R-CDDU/julia_code";

function readin_file(folderPath, n)
    data_m = Dict();
    for i4 in fileList
        # read the data from the csv files
        myPath = folderPath * i4 * ".csv";
        data_m[i4[2:length(i4)]] = readdlm(myPath,',');
    end
    # convert the raw data to the data structure
    T = 9;

    # market data
    D = data_m["Market"][1,1:T];
    h = data_m["Market"][2,1:T];
    s = data_m["Market"][3,1:T];

    if mod(n,3) == 2
        n1 = Int64(floor(n/3)+1);
        n2 = Int64(floor(n/3)+1);
        n3 = n - n1 - n2;
    elseif mod(n,3) == 1
        n1 = Int64(floor(n/3)+1);
        n2 = Int64(floor(n/3));
        n3 = n - n1 - n2;
    else
        n1 = Int64(n/3);
        n2 = Int64(n/3);
        n3 = Int64(n/3);
    end

    # resource data
    cc = vcat(data_m["Resource1"][1,1:n1], data_m["Resource2"][1,1:n2], data_m["Resource3"][1,1:n3]);
    RD = vcat(data_m["Resource1"][2,1:n1], data_m["Resource2"][2,1:n2], data_m["Resource3"][2,1:n3]);
    RU = vcat(data_m["Resource1"][3,1:n1], data_m["Resource2"][3,1:n2], data_m["Resource3"][3,1:n3]);
    pmax = vcat(data_m["Resource1"][4,1:n1], data_m["Resource2"][4,1:n2], data_m["Resource3"][4,1:n3]);
    pmin = vcat(data_m["Resource1"][5,1:n1], data_m["Resource2"][5,1:n2], data_m["Resource3"][5,1:n3]);
    upmin = vcat(data_m["Resource1"][6,1:n1], data_m["Resource2"][6,1:n2], data_m["Resource3"][6,1:n3]);
    downmin = vcat(data_m["Resource1"][7,1:n1], data_m["Resource2"][7,1:n2], data_m["Resource3"][7,1:n3]);
    alpha = vcat(data_m["Resource1"][8,1:n1], data_m["Resource2"][8,1:n2], data_m["Resource3"][8,1:n3]);
    beta = vcat(data_m["Resource1"][9,1:n1], data_m["Resource2"][9,1:n2], data_m["Resource3"][9,1:n3]);

    mData = DR_data(T, D, h, s, cc, RD, RU, pmax, pmin, upmin, downmin, alpha, beta);
    return mData;
end

data = Dict();
for i1 in firstsub
    i1l = length(i1);
    data[i1[i1l-1:i1l]] = Dict();
    for i2 in secondsub
        i2l = length(i2);
        data[i1[i1l-1:i1l]][i2[i2l-1:i2l]] = Dict();
        for i3ind in 1:length(lastsub)
            i3 = lastsub[i3ind];
            i3n = lastn[i3ind];
            try
                mData = readin_file(rootPath* i1 * i2 * i3, i3n);
                data[i1[i1l-1:i1l]][i2[i2l-1:i2l]][i3[2:length(i3)]] = mData;
            catch
                a = 1;
            end
        end
    end
end

save("testData.jld","mData",data);
