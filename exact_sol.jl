# obtain an exact solution
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

data_raw = load("testData.jld");
data = data_raw["mData"];
mData = deepcopy(data["30"]["30"]["D1N1200"]);
clusterList, clusterInfo, clusterRev = genClusters(1200,24)

Gamma = Dict();
for t in 1:mData.T
    Gamma[t] = 0.07;
end

n = 1200;
mData = combineGen(mData);

# compare mccormick with bilinear for the gap
mcp_r = mccormick_opt_relax(mData, n, clusterList, clusterRev, Gamma);
optimize!(mcp_r);

mcp = mccormick_opt(mData, n, clusterList, clusterRev, Gamma);
optimize!(mcp);

bps = bilinear_simple(mData, n, clusterList, clusterRev, Gamma, 3600);
optimize!(bps);
# MIP
mp_mip = mip_form(mData, n, clusterList, clusterRev, Gamma, 3600);
optimize!(mp_mip);
# ADA
vbest_ada, xbest_ada = main_ada(mData, n, clusterList, clusterRev, Gamma);
# CGA
vbest_cga, xbest_cga = main_cga(mData, n, clusterList, clusterRev, Gamma);
