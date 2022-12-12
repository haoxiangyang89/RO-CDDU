# deterministic tests
using DelimitedFiles, JLD, HDF5;
using JuMP, Gurobi, CSV;
const GUROBI_ENV = Gurobi.Env();

include("def.jl");
include("auxiliary.jl");
include("det.jl");

data_raw = load("testData.jld");
data = data_raw["mData"];
n = 800;
cl_no = 20;

for tolerance_level in ["10","30"]
    for demand_level in ["30","80"]
        mData = data[tolerance_level][demand_level]["D1N$(n)"];
        clusterList, clusterInfo, clusterRev = genClusters(n,cl_no);
        mp2 = det_opt_2(mData, n, clusterList, clusterRev);
        startTime = time();
        optimize!(mp2);
        time_det = time() - startTime;
        x_det = Dict();
        for i in 1:n
            for t in 1:mData.T
                x_det[i,t] = value(mp2[:x][i,t]);
            end
        end
        v_det = objective_value(mp2);
        save("det_time_results_$(tolerance_level)_$(demand_level)_$(n).jld", "time", time_det,
                                          "obj", v_det,
                                          "x_best", x_det);

    end
end

#mp = det_opt(mData, 1200);
