# cga test
using DelimitedFiles, JLD, HDF5;
using JuMP, Gurobi, CSV;
const GUROBI_ENV = Gurobi.Env();

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

# nList = [5, 50, 200, 400, 800, 1200];
# clusterNo = [1, 5, 20, 20, 20, 20];

tolerance_level = ARGS[1];
demand_level = ARGS[2];

n = parse(Int64, ARGS[3]);
cl_no = parse(Int64, ARGS[4]);

mData = deepcopy(data[tolerance_level][demand_level]["D1N$(n)"]);
clusterList, clusterInfo, clusterRev = genClusters(n, cl_no);

# for shortage tolerance 30, demand level 30, Gamma = 0.07
gamma = parse(Float64, ARGS[5]);
Gamma = Dict();
for t in 1:mData.T
    Gamma[t] = gamma;
end

# CGA
time_cga, vbest_cga, xbest_cga, πbest_cga = main_cga(mData, n, clusterList, clusterRev, Gamma, 1, 36000);
save("cga_time_results_$(tolerance_level)_$(demand_level)_$(n)_$(gamma*100).jld", "time", time_cga,
                                  "obj", vbest_cga,
                                  "x_best", xbest_cga,
                                  "pi_best", πbest_cga);

println("======================================= CGA Case $(tolerance_level)_$(demand_level)_$(n)_$(gamma) Finished =======================================");
println("================= Time $(time_cga), Best Value $(vbest_cga) =================");
