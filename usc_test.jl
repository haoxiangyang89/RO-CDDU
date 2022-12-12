# test the algorithm against the USC data
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
data_raw = load("testData_USC.jld");
data = data_raw["mData"];
mData = deepcopy(data["1400"]);

n = 115;
cl_no = 115;
clusterList, clusterInfo, clusterRev = genClusters(n, cl_no);

gamma = 0.01;
Gamma = Dict();
for t in 1:mData.T
    Gamma[t] = gamma;
end

time_limit = 36000;
# test the mip result
# mip_status, time_mip, mip_lb, mip_ub, xDict_mip, πDict_mip = solve_mip(mData, n, clusterList, clusterRev, Gamma, time_limit);
# bp_status, time_bp, bp_lb, bp_ub, xDict_bp, πDict_bp = solve_bilinear(mData, n, clusterList, clusterRev, Gamma, 360);
# test the cga result
# time_cga, vbest_cga, xbest_cga, πbest_cga = main_cga(mData, n, clusterList, clusterRev, Gamma, 1, 36000);

gamma_List = 0.01*[1,2,3,4,5,6,7,8,9,10];
for gamma in gamma_List
    Gamma = Dict();
    for t in 1:mData.T
        Gamma[t] = gamma;
    end
    πhat = Dict();
    πhat[1] = Dict();
    πhat[2] = Dict();
    πhat[4] = Dict();
    for t in 1:mData.T
        πhat[4][t,1] = 0.0;
        πhat[4][t,2] = 0.0;
        for i in 1:n
            πhat[1][i,t] = mData.h[t] - mData.cc[i];
            πhat[2][i,t] = -mData.s[t] - mData.cc[i];
        end
    end

    # initialize another π
    πhat2 = Dict();
    πhat2[1] = Dict();
    πhat2[2] = Dict();
    πhat2[4] = Dict();
    for t in 1:mData.T
        πhat2[4][t,1] = maximum([mData.h[t] - mData.cc[i] for i in 1:n]);
        πhat2[4][t,2] = maximum([mData.s[t] + mData.cc[i] for i in 1:n]);
        for i in 1:n
            πhat2[1][i,t] = 0.0;
            πhat2[2][i,t] = 0.0;
        end
    end

    πList = [πhat, πhat2];

    # initialize the big-M
    M1,M2 = getM(n, mData, Gamma);

    # initialize the xprob
    xprob = cga_ini_xprob(mData, n, clusterList, clusterRev, Gamma);
    xhat, vbest, μhat = solve_x(mData, xprob, πList, n, clusterList, clusterRev, Gamma, M1, M2);
    xbest = deepcopy(xhat);
    πbest = Dict();
    μbest = Dict();

    stopBool = false;
    start_time = time();
    while !(stopBool)
        πhat = solve_pi(mData, xhat, n, clusterList, clusterRev, Gamma);
        push!(πList, πhat);
        μhat_new = [];
        if mode == 1
            xhat_new, vhat_new, μhat_new = solve_x(mData, xprob, πList, n, clusterList, clusterRev, Gamma, M1, M2, 0, 10800, xhat, μhat);
        else
            xhat_new, vhat_new = solve_x_new(mData, xprob, πList, n, clusterList, clusterRev, Gamma);
        end

        if sum(sum(abs(xhat[i,t] - xhat_new[i,t]) for t in 1:mData.T) for i in 1:n) <= 0.1
            stopBool = true;
        else
            xhat = deepcopy(xhat_new);
            elapsed_time = time() - start_time;
            if elapsed_time > time_limit
                stopBool = true;
            end
            μhat = deepcopy(μhat_new);
        end
        if vhat_new < vbest
            vbest = vhat_new;
            xbest = xhat_new;
            μbest = deepcopy(μhat_new);
        end
    end

    πbest = solve_pi(mData, xbest, n, clusterList, clusterRev, Gamma);
    elapsed_time = time() - start_time;

    save("cga_time_results_usc_$(gamma*100).jld", "time", elapsed_time,
                                      "obj", vbest,
                                      "x_best", xbest,
                                      "pi_best", πbest);
end

data_result = Dict();
for i in 1:10
    # read in the result data
    data_result = load("cga_time_results_usc_$(i).0.jld");
    # output the x in csv format
    mat_x = zeros(n,mData.T);
    for ni in 1:n
        for t in 1:mData.T
            mat_x[ni,t] = data_result["x_best"][ni,t];
        end
    end
    writedlm("usc_x_$(i).csv",mat_x,",")
end
