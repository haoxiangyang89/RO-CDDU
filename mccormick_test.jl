# McCormick test
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

# McCormick relaxation
mcr_status, time_mcr, mcr_lb, mcr_ub, xDict_mcr, πDict_mcr = solve_mccormick(mData, n, clusterList, clusterRev, Gamma, 1, 36000);
save("mcr_time_results_$(tolerance_level)_$(demand_level)_$(n)_$(gamma*100).jld", "status", mcr_status,
                                  "time", time_mcr,
                                  "lb", mcr_lb,
                                  "ub", mcr_ub,
                                  "x_best", xDict_mcr,
                                  "pi_best", πDict_mcr);

println("======================================= McCormick LP Relaxation Case $(tolerance_level)_$(demand_level)_$(n)_$(gamma) Finished =======================================");
println("================= Status $(mcr_status), Time $(time_mcr), Best Value $(mcr_ub) =================");


# McCormick relaxation
mcp_status, time_mcp, mcp_lb, mcp_ub, xDict_mcp, πDict_mcp = solve_mccormick(mData, n, clusterList, clusterRev, Gamma, 2, 36000);
save("mcp_time_results_$(tolerance_level)_$(demand_level)_$(n)_$(gamma*100).jld", "status", mcp_status,
                                "time", time_mcp,
                                "lb", mcp_lb,
                                "ub", mcp_ub,
                                "x_best", xDict_mcp,
                                "pi_best", πDict_mcp);

println("======================================= McCormick MIP Relaxation Case $(tolerance_level)_$(demand_level)_$(n)_$(gamma) Finished =======================================");
println("================= Status $(mcp_status), Time $(time_mcp), Best Value $(mcp_ub) =================");
