# solution tests: cga vs. ada for special case in Appendix 3
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

tolerance_level = "30";
demand_level = "30";

m = 1200;
mData = deepcopy(data[tolerance_level][demand_level]["D1N$(m)"]);

n = 5;
cl_no = 1;
clusterList, clusterInfo, clusterRev = genClusters(n, cl_no);
for t in 1:mData.T
   mData.D[t] = mData.D[t]/(m/n);
end

gamma = 0.07;
Gamma = Dict();
for t in 1:mData.T
    Gamma[t] = gamma;
end

time_ada, vbest_ada, xbest_ada, πbest_ada = main_ada(mData, n, clusterList, clusterRev, Gamma, 36000);
time_cga, vbest_cga, xbest_cga, πbest_cga = main_cga(mData, n, clusterList, clusterRev, Gamma, 36000);
mip_status, time_mip, mip_lb, mip_ub, xDict_mip, πDict_mip = solve_mip(mData, n, clusterList, clusterRev, Gamma, 18000);
