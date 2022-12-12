# generate random scenarios

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

scenTot = Dict();
for tolerance_level in ["10", "30"]
    scenTot[tolerance_level] = Dict();
    for demand_level in ["30", "80"]
        # read the test data for n = 800
        mData = deepcopy(data[tolerance_level][demand_level]["D1N800"]);
        scenTot[tolerance_level][demand_level] = scenGen(mData, 800, 5000);
    end
end

save("scenData.jld","scenTot",scenTot);
