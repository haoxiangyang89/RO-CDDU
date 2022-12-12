# load gamma results and print out results
using DelimitedFiles, JLD, HDF5;
using JuMP, Gurobi, CSV;

include("def.jl");
include("auxiliary.jl");
include("mccormick.jl");
include("bilinear.jl");
include("ada.jl");
include("cga.jl");
include("mip_form.jl");

# load the data
data_raw = load("testData.jld");
data = data_raw["mData"];

nList = [5, 50, 200, 400, 800, 1200];
clusterNo = [1, 5, 20, 20, 20, 20];

objCGA = Dict();
timeCGA = Dict();
#n = 20;
n = 50;
for tolerance_level in ["30","10"]
    for demand_level in ["30","80"]
        objCGA[tolerance_level,demand_level] = Dict();
        timeCGA[tolerance_level,demand_level] = Dict();
        for gamma in 1:10
            data_cga = load("./test_results_jlds/Gamma_$(n)/cga_time_results_$(tolerance_level)_$(demand_level)_$(n)_$(gamma).0.jld");
            objCGA[tolerance_level,demand_level][gamma] = data_cga["obj"];
            timeCGA[tolerance_level,demand_level][gamma] = data_cga["time"];
        end
    end
end

for tolerance_level in ["30","10"]
    for demand_level in ["30","80"]
        println("------ CGA Case $(tolerance_level), $(demand_level) -----")
        for gamma in 1:10
            print("& $(-round(objCGA[tolerance_level,demand_level][gamma],digits = 2)) ");
        end
        println("\\\\");
        for gamma in 1:10
            print("& $(round(timeCGA[tolerance_level,demand_level][gamma],digits = 1)) ");
        end
        println("\\\\");
    end
end

objADA = Dict();
timeADA = Dict();
for tolerance_level in ["30","10"]
    for demand_level in ["30","80"]
        objADA[tolerance_level,demand_level] = Dict();
        timeADA[tolerance_level,demand_level] = Dict();
        for gamma in 1:10
            data_ada = load("./test_results_jlds/Gamma_$(n)/ada_time_results_$(tolerance_level)_$(demand_level)_$(n)_$(gamma).0.jld");
            objADA[tolerance_level,demand_level][gamma] = data_ada["obj"];
            timeADA[tolerance_level,demand_level][gamma] = data_ada["time"];
        end
    end
end

for tolerance_level in ["30","10"]
    for demand_level in ["30","80"]
        println("------ ADA Case $(tolerance_level), $(demand_level) -----")
        for gamma in 1:10
            print("& $(-round(objADA[tolerance_level,demand_level][gamma],digits = 2)) ");
        end
        println("\\\\");
        for gamma in 1:10
            print("& $(round(timeADA[tolerance_level,demand_level][gamma],digits = 1)) ");
        end
        println("\\\\");
    end
end

objMIP = Dict();
timeMIP = Dict();
lbMIP = Dict();
for tolerance_level in ["30","10"]
    for demand_level in ["30","80"]
        objMIP[tolerance_level,demand_level] = Dict();
        timeMIP[tolerance_level,demand_level] = Dict();
        lbMIP[tolerance_level,demand_level] = Dict();
        for gamma in 1:10
            data_mip = load("./test_results_jlds/Gamma_$(n)/mip_time_results_$(tolerance_level)_$(demand_level)_$(n)_$(gamma).0.jld");
            objMIP[tolerance_level,demand_level][gamma] = data_mip["ub"];
            lbMIP[tolerance_level,demand_level][gamma] = data_mip["lb"];
            timeMIP[tolerance_level,demand_level][gamma] = data_mip["time"];
        end
    end
end

for tolerance_level in ["30","10"]
    for demand_level in ["30","80"]
        println("------ MIP Case $(tolerance_level), $(demand_level) -----")
        for gamma in 1:10
            print("& $(-round(objMIP[tolerance_level,demand_level][gamma],digits = 2)) ");
        end
        println("\\\\");
        for gamma in 1:10
            print("& $(round(timeMIP[tolerance_level,demand_level][gamma],digits = 1)) ");
        end
        println("\\\\");
    end
end
