# simulation analysis
using DelimitedFiles, JLD, HDF5, Statistics;
using JuMP, Gurobi, CSV;
const GUROBI_ENV = Gurobi.Env();

include("def.jl");
include("det.jl");
include("auxiliary.jl");
include("mccormick.jl");
include("bilinear.jl");
include("ada.jl");
include("cga.jl");
include("mip_form.jl");

# load the data
data_raw = load("testData.jld");
data = data_raw["mData"];
scen_raw = load("scenData.jld");
scenTot = scen_raw["scenTot"];

# nList = [5, 50, 200, 400, 800, 1200];
# clusterNo = [1, 5, 20, 20, 20, 20];

n = 800;
solution_analysis_data = Dict();

for tolerance_level in ["10","30"]
    for demand_level in ["30", "80"]
        gen_time = Dict();
        simu_data = Dict();
        gen_proportion = Dict();

        mData = deepcopy(data[tolerance_level][demand_level]["D1N$(n)"]);
        scenData = deepcopy(scenTot[tolerance_level][demand_level]);

        # obtain the solution from generating the solution file
        # first load the solution file
        result_data = load("./test_results_jlds/Gamma_800/det_time_results_$(tolerance_level)_$(demand_level)_800.jld");
        xbest = result_data["x_best"];
        # obtain the simulation result
        randList = zeros(5000);
        for ω in 1:5000
            randList[ω] = eval_opt(mData, n, xbest, scenData[ω]);
        end
        simu_data[0] = randList;
        gen_time[0] = [sum(xbest[i,t] for i in 1:n) for t in 1:mData.T];

        # obtain the type a/b/c generation percentage
        totalA = sum(sum(xbest[i,t] for i in 1:267) for t in 1:mData.T);
        totalB = sum(sum(xbest[i,t] for i in 268:524) for t in 1:mData.T);
        totalC = sum(sum(xbest[i,t] for i in 525:800) for t in 1:mData.T);
        percentA = totalA/(totalA + totalB + totalC);
        percentB = totalB/(totalA + totalB + totalC);
        percentC = totalC/(totalA + totalB + totalC);
        gen_proportion[0] = [[totalA, totalB, totalC], [percentA, percentB, percentC]];

        println("======================================= Simulation Case $(tolerance_level)_$(demand_level)_$(n)_0 Finished =======================================");
        println("================= Average Cost $(mean(randList)), Robust Cost $(result_data["obj"]) =================");
        println("================= Type A Generation Amount: $(totalA) ($(round(percentA*100,digits = 2))%), Type B Generation Amount: $(totalB) ($(round(percentB*100,digits = 2))%), Type C Generation Amount: $(totalC) ($(round(percentC*100,digits = 2))%)=================");
        println("");

        for gInd in 1:10
            result_data = load("./test_results_jlds/Gamma_800/cga_time_results_$(tolerance_level)_$(demand_level)_800_$(gInd).0.jld");
            xbest = result_data["x_best"];

            # obtain the simulation result
            randList = zeros(5000);
            for ω in 1:5000
                randList[ω] = eval_opt(mData, n, xbest, scenData[ω]);
            end
            simu_data[gInd] = randList;
            gen_time[gInd] = [sum(xbest[i,t] for i in 1:n) for t in 1:mData.T];

            # obtain the type a/b/c generation percentage
            totalA = sum(sum(xbest[i,t] for i in 1:267) for t in 1:mData.T);
            totalB = sum(sum(xbest[i,t] for i in 268:524) for t in 1:mData.T);
            totalC = sum(sum(xbest[i,t] for i in 525:800) for t in 1:mData.T);
            percentA = totalA/(totalA + totalB + totalC);
            percentB = totalB/(totalA + totalB + totalC);
            percentC = totalC/(totalA + totalB + totalC);
            gen_proportion[gInd] = [[totalA, totalB, totalC], [percentA, percentB, percentC]];

            println("======================================= Simulation Case $(tolerance_level)_$(demand_level)_$(n)_$(gInd*0.01) Finished =======================================");
            println("================= Average Cost $(mean(randList)), Robust Cost $(result_data["obj"]) =================");
            println("================= Type A Generation Amount: $(totalA) ($(round(percentA*100,digits = 2))%), Type B Generation Amount: $(totalB) ($(round(percentB*100,digits = 2))%), Type C Generation Amount: $(totalC) ($(round(percentC*100,digits = 2))%)=================");
            println("");
        end

        # record the data
        solution_analysis_data[tolerance_level, demand_level] = [simu_data, gen_time, gen_proportion];
    end
end
save("sol_analysis_data.jld", "sol_data", solution_analysis_data);

# print out the csv files for ABC proportions
for tolerance_level in ["10","30"]
    for demand_level in ["30", "80"]
        print_ABC = zeros(11,4);
        for gInd in 0:10
            print_ABC[gInd + 1, 1] = gInd;
            print_ABC[gInd + 1, 2] = solution_analysis_data[tolerance_level, demand_level][3][gInd][2][1];
            print_ABC[gInd + 1, 3] = solution_analysis_data[tolerance_level, demand_level][3][gInd][2][2];
            print_ABC[gInd + 1, 4] = solution_analysis_data[tolerance_level, demand_level][3][gInd][2][3];
        end
        open("proportion_$(tolerance_level)_$(demand_level).csv", "w") do io
             writedlm(io, print_ABC, ',');
        end
    end
end

# print out the csv files for demand response amount
for tolerance_level in ["10","30"]
    for demand_level in ["30", "80"]
        nom_demand = data[tolerance_level][demand_level]["D1N$(n)"].D;
        print_DR = zeros(9,6);
        for tInd in 1:9
            print_DR[tInd, 1] = tInd;
            print_DR[tInd, 2] = nom_demand[tInd];
            print_DR[tInd, 3] = solution_analysis_data[tolerance_level, demand_level][2][0][tInd];
            print_DR[tInd, 4] = solution_analysis_data[tolerance_level, demand_level][2][1][tInd];
            print_DR[tInd, 5] = solution_analysis_data[tolerance_level, demand_level][2][5][tInd];
            print_DR[tInd, 6] = solution_analysis_data[tolerance_level, demand_level][2][10][tInd];
        end
        open("DRamount_$(tolerance_level)_$(demand_level).csv", "w") do io
             writedlm(io, print_DR, ',');
        end
    end
end

# print out the csv files for simulated expectation results
for tolerance_level in ["10","30"]
    for demand_level in ["30", "80"]
        print_simu = zeros(11, 4);
        for gInd in 0:10
            meanV = -mean(solution_analysis_data[tolerance_level, demand_level][1][gInd]);
            stdV = std(solution_analysis_data[tolerance_level, demand_level][1][gInd]);
            print_simu[gInd + 1, 1] = gInd;
            print_simu[gInd + 1, 2] = meanV;
            print_simu[gInd + 1, 3] = meanV + 1.96 * stdV;
            print_simu[gInd + 1, 4] = meanV - 1.96 * stdV;
        end
        open("simu_$(tolerance_level)_$(demand_level).csv", "w") do io
             writedlm(io, print_simu, ',');
        end
    end
end
