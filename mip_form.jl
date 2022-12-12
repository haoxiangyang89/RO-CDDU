# mip full formulation
function mip_form(mData, n, clusterList, clusterRev, Gamma, time_limit = 3600)
    # obtain the big-M
    M1, M2 = getM(n, mData, Gamma);

    # enumerate the πList
    πDict = Dict();
    πDict[1] = Dict();
    πDict[2] = Dict();
    for t in 1:mData.T
        πDict[1][t] = [];
        πDict[2][t] = [];
    end
    xmax = sum(mData.pmax[i] for i in 1:n);

    for t in 1:mData.T
        hcList = sort(unique([mData.h[t] - mData.cc[i] for i in 1:n]));
        pushfirst!(hcList,0.0);
        for hcInd in 1:length(hcList)
            π1hat = Dict(1 => [max(mData.h[t] - mData.cc[i] - hcList[hcInd], 0) for i in 1:n],
                         4 => hcList[hcInd]);
            push!(πDict[1][t], π1hat)
        end
    end

    for t in 1:mData.T
        scList = sort(unique([mData.s[t] + mData.cc[i] for i in 1:n]));
        pushfirst!(scList,0.0);
        for scInd in 1:length(scList)
            π2hat = Dict(2 => [min(-mData.s[t] - mData.cc[i] + scList[scInd], 0) for i in 1:n],
                         4 => scList[scInd]);
            push!(πDict[2][t], π2hat)
        end
    end

    mp = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(GUROBI_ENV),
                                         "TimeLimit" => time_limit,
                                         "Threads" => 30));
    # set up the variables:
    #   x: amount of demand response
    #   u: commitment indicator
    #   w,v: increase/decrease indicator
    @variable(mp, x[i in 1:n, t in 1:mData.T]);
    @variable(mp, zplus[t in 1:mData.T] >= 0);
    @variable(mp, zminus[t in 1:mData.T] >= 0);
    @variable(mp, u[i in clusterList, t in 1:mData.T], Bin);
    @variable(mp, w[i in clusterList, t in 2:mData.T], Bin);
    @variable(mp, v[i in clusterList, t in 2:mData.T], Bin);

    # dual variables
    @variable(mp, π1[1:n, 1:mData.T] >= 0);
    @variable(mp, π2[1:n, 1:mData.T] <= 0);
    @variable(mp, π4[1:mData.T, 1:2] >= 0);

    @variable(mp, V[t in 1:mData.T]);
    @constraint(mp, V1[t in 1:mData.T], V[t] >= sum((mData.h[t] - mData.cc[i]) * x[i,t] for i in 1:n) -
                                                mData.h[t] * mData.D[t] + zplus[t]);
    @constraint(mp, V2[t in 1:mData.T], V[t] >= sum((-mData.s[t] - mData.cc[i]) * x[i,t] for i in 1:n) +
                                                mData.s[t] * mData.D[t] + zminus[t]);
    # commitment bounds
    @constraint(mp, commitment_up[i in 1:n, t in 1:mData.T], x[i,t] <= mData.pmax[i] * u[clusterRev[i],t]);
    @constraint(mp, commitment_down[i in 1:n, t in 1:mData.T], x[i,t] >= mData.pmin[i] * u[clusterRev[i],t]);

    # ramping constraints
    @constraint(mp, ramp_up[i in 1:n, t in 1:(mData.T - 1)], x[i,t + 1] - x[i,t] <= mData.RU[i] * u[clusterRev[i],t + 1]);
    @constraint(mp, ramp_down[i in 1:n, t in 1:(mData.T - 1)], x[i,t + 1] - x[i,t] >= -mData.RD[i] * u[clusterRev[i],t]);

    # cooling down constraints
    @constraint(mp, w1[i in 1:n, t in 2:mData.T], x[i,t] - x[i,t - 1] <= mData.RU[i] * w[clusterRev[i],t]);
    @constraint(mp, w2[i in 1:n, t in 2:mData.T], x[i,t] - x[i,t - 1] >= -mData.RD[i] * v[clusterRev[i],t]);
    @constraint(mp, w3[i in 1:n, t in 2:mData.T, τ in t:min(t + Int64(mData.upmin[i]) - 1,mData.T)], x[i,τ] - x[i,τ - 1] >= -mData.RD[i] * (1 - w[clusterRev[i],t]));
    @constraint(mp, w4[i in 1:n, t in 2:mData.T, τ in t:min(t + Int64(mData.downmin[i]) - 1,mData.T)], x[i,τ] - x[i,τ - 1] <= mData.RU[i] * (1 - v[clusterRev[i],t]));

    # tightening constraint
    @constraint(mp, tighten1[i in clusterList, t in 2:mData.T], w[i,t] <= u[i,t]);
    @constraint(mp, tighten2[i in clusterList, t in 2:mData.T], v[i,t] <= u[i,t - 1]);
    @constraint(mp, tighten3[i in clusterList, t in 2:mData.T], w[i,t] + v[i,t] <= 1);

    # add all possible π's
    @variable(mp, μ[t in 1:mData.T, l in 1:2, cs in 1:length(πDict[l][t])], Bin);
    @constraint(mp, μCons[t in 1:mData.T, l in 1:2], sum(μ[t, l, cs] for cs in 1:length(πDict[l][t])) == 1);

    # add the current columns
    @constraint(mp, zplusCons[t in 1:mData.T, cs in 1:length(πDict[1][t])], mp[:zplus][t] >= sum(πDict[1][t][cs][1][i] * mData.beta[i] * mp[:x][i,t] for i in 1:n) +
                        πDict[1][t][cs][4] * Gamma[t] * sum(mData.pmax[i] for i in 1:n) - M1[t] * (1 - μ[t,1,cs]));
    @constraint(mp, zminusCons[t in 1:mData.T, cs in 1:length(πDict[2][t])], mp[:zminus][t] >= sum(πDict[2][t][cs][2][i] * mData.alpha[i] * mp[:x][i,t] for i in 1:n) +
                        πDict[2][t][cs][4] * Gamma[t] * sum(mData.pmax[i] for i in 1:n) - M2[t] * (1 - μ[t,2,cs]));

    @objective(mp, Min, sum(V[t] for t in 1:mData.T));
    return mp;
end

function solve_mip(mData, n, clusterList, clusterRev, Gamma, time_limit)
    # create/solve a bilinear program
    mip_mp = mip_form(mData, n, clusterList, clusterRev, Gamma, time_limit);
    start_time = time();
    optimize!(mip_mp);
    elapsed_time = time() - start_time;

    # return status, solution, and optimal value
    mip_mp_status = termination_status(mip_mp);
    mip_mp_lb = objective_bound(mip_mp);
    mip_mp_ub = objective_value(mip_mp);
    x_val = value.(mip_mp[:x]);
    π1_val = value.(mip_mp[:π1]);
    π2_val = value.(mip_mp[:π2]);
    π4_val = value.(mip_mp[:π4]);
    xDict = Dict();
    πDict = Dict(1 => Dict(),
                 2 => Dict(),
                 4 => Dict());

    for t in 1:mData.T
        for i in 1:n
            xDict[i,t] = x_val[i,t];
            πDict[1][i,t] = π1_val[i,t];
            πDict[2][i,t] = π2_val[i,t];
        end
        πDict[4][t,1] = π4_val[t,1];
        πDict[4][t,2] = π4_val[t,2];
    end
    return mip_mp_status, elapsed_time, mip_mp_lb, mip_mp_ub, xDict, πDict;
end
