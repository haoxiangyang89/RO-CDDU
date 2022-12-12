# test with different n
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

n = 2;

Gamma = Dict();
for t in 1:mData.T
    Gamma[t] = 0.07;
end

clusterList, clusterInfo, clusterRev = genClusters(2,1)
for t in 1:mData.T
    mData.D[t] = mData.D[t]/(1200/n);
end

# MIP
mp_mip = mip_form(mData, n, clusterList, clusterRev, Gamma, 3600);
optimize!(mp_mip);
# ADA
vbest_ada, xbest_ada = main_ada(mData, n, clusterList, clusterRev, Gamma);
# CGA
vbest_cga, xbest_cga = main_cga(mData, n, clusterList, clusterRev, Gamma);
