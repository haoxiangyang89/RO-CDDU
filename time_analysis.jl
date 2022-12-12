# load test results and print out results
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

nList = [5, 50, 200, 400, 800, 1200];
clusterNo = [1, 5, 20, 20, 20, 20];

for nInd in 1:length(nList)
    n = nList[nInd];

    data_mcr = load("./test_results_jlds/mcr_time_results_30_30_$(n)_5.0.jld");
    println("Case $(n), Status $(data_mcr["status"]), Time $(data_mcr["time"]), Obj $(data_mcr["ub"])");
end

dataMCP = Dict();
for nInd in 1:length(nList)
    n = nList[nInd];

    data_mcp = load("./test_results_jlds/mcp_time_results_30_30_$(n)_5.0.jld");
    println("Case $(n), Status $(data_mcp["status"]), Time $(data_mcp["time"]), Obj $(data_mcp["ub"])");
    dataMCP[n] = data_mcp["ub"];
end

dataADA = Dict();
for nInd in 1:length(nList)
    n = nList[nInd];

    data_ada = load("./test_results_jlds/ada_time_results_30_30_$(n)_5.0.jld");
    println("Case $(n), Time $(data_ada["time"]), Obj $(data_ada["obj"])");
    dataADA[n] = data_ada["obj"];
end

dataCGA = Dict();
for nInd in 1:length(nList)
    n = nList[nInd];

    data_cga = load("./test_results_jlds/cga_time_results_30_30_$(n)_5.0.jld");
    println("Case $(n), Time $(data_cga["time"]), Obj $(data_cga["obj"])");
    dataCGA[n] = data_cga["obj"];
end

dataMIP = Dict();
for nInd in 1:length(nList)
    n = nList[nInd];

    data_mip = load("./test_results_jlds/mip_time_results_30_30_$(n).jld");
    println("Case $(n), Status $(data_mip["status"]), Time $(data_mip["time"]), Obj $(data_mip["ub"]), LB $(data_mip["lb"])");
    dataMIP[n] = (data_mip["ub"], data_mip["lb"]);
end

for n in nList
    print("& $(round(-dataMIP[n][1]), digits = 2) ");
end
println("\\\\");
for n in nList
    print("& $(round((dataMIP[n][1] - dataMIP[n][2])/(-dataMIP[n][2]),digits = 2)) ");
end
println("\\\\");
