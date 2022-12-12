# mccormick relaxation
function mccormick_opt_relax(mData, n, clusterList, clusterRev, Gamma, time_limit = 3600)
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
    @variable(mp, 0 <= u[i in clusterList, t in 1:mData.T] <= 1);
    @variable(mp, 0 <= w[i in clusterList, t in 2:mData.T] <= 1);
    @variable(mp, 0 <= v[i in clusterList, t in 2:mData.T] <= 1);

    # dual variables
    @variable(mp, π1[1:n, 1:mData.T] >= 0);
    @variable(mp, π2[1:n, 1:mData.T] <= 0);
    @variable(mp, π4[1:mData.T, 1:2] >= 0);

    # π*u variables
    @variable(mp, pu1[1:n, 1:mData.T] >= 0);
    @variable(mp, pu2[1:n, 1:mData.T] <= 0);

    # mccormick relaxation terms
    @variable(mp, q1[1:n, 1:mData.T]);
    @variable(mp, q2[1:n, 1:mData.T]);

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

    # bilinear reformulation
    @constraint(mp, zplusCons[t in 1:mData.T], zplus[t] >= sum(q1[i,t] * mData.beta[i] for i in 1:n) +
                        π4[t,1] * Gamma[t] * sum(mData.pmax[i] for i in 1:n));
    @constraint(mp, zminusCons[t in 1:mData.T], zminus[t] >= sum(q2[i,t] * mData.alpha[i] for i in 1:n) +
                        π4[t,2] * Gamma[t] * sum(mData.pmax[i] for i in 1:n));

    # dual feasibility constraints
    @constraint(mp, dualCons1_plus[i in 1:n, t in 1:mData.T], mData.h[t] - mData.cc[i] - π1[i,t] <= π4[t,1]);
    @constraint(mp, dualCons2_plus[i in 1:n, t in 1:mData.T], -(mData.h[t] - mData.cc[i]) + π1[i,t] <= π4[t,1]);
    @constraint(mp, dualCons1_minus[i in 1:n, t in 1:mData.T], (-mData.s[t] - mData.cc[i]) - π2[i,t] <= π4[t,2]);
    @constraint(mp, dualCons2_minus[i in 1:n, t in 1:mData.T], mData.s[t] + mData.cc[i] + π2[i,t] <= π4[t,2]);

    # mccormick relaxation
    πbar1 = Dict();
    πubar2 = Dict();
    for i in 1:n
        for t in 1:mData.T
            πbar1[i,t] = mData.h[t] - mData.cc[i];
            πubar2[i,t] = -mData.s[t] - mData.cc[i];
        end
    end

    # π*u reformulation
    @constraint(mp, pu1Cons1[i in 1:n, t in 1:mData.T], pu1[i,t] <= π1[i,t]);
    @constraint(mp, pu1Cons2[i in 1:n, t in 1:mData.T], pu1[i,t] <= (mData.h[t] - mData.cc[i]) * u[clusterRev[i],t]);
    @constraint(mp, pu1Cons3[i in 1:n, t in 1:mData.T], pu1[i,t] >= π1[i,t] - (mData.h[t] - mData.cc[i]) * (1 - u[clusterRev[i],t]));

    @constraint(mp, pu2Cons1[i in 1:n, t in 1:mData.T], pu2[i,t] >= π2[i,t]);
    @constraint(mp, pu2Cons2[i in 1:n, t in 1:mData.T], pu2[i,t] >= (-mData.s[t] - mData.cc[i]) * u[clusterRev[i],t]);
    @constraint(mp, pu2Cons3[i in 1:n, t in 1:mData.T], pu2[i,t] <= π2[i,t] - (-mData.s[t] - mData.cc[i]) * (1 - u[clusterRev[i],t]));

    # McCormick relaxationn
    @constraint(mp, q1Cons1[i in 1:n, t in 1:mData.T], q1[i,t] >= pu1[i,t] * mData.pmin[i]);
    @constraint(mp, q1Cons2[i in 1:n, t in 1:mData.T], q1[i,t] >= pu1[i,t] * mData.pmax[i] + x[i,t] * (mData.h[t] - mData.cc[i]) -
                                                                    (mData.h[t] - mData.cc[i]) * mData.pmax[i] * u[clusterRev[i],t]);
    @constraint(mp, q1Cons3[i in 1:n, t in 1:mData.T], q1[i,t] <= pu1[i,t] * mData.pmin[i] + x[i,t] * (mData.h[t] - mData.cc[i]) -
                                                                    (mData.h[t] - mData.cc[i]) * mData.pmin[i] * u[clusterRev[i],t]);
    @constraint(mp, q1Cons4[i in 1:n, t in 1:mData.T], q1[i,t] <= pu1[i,t] * mData.pmax[i]);

    @constraint(mp, q2Cons1[i in 1:n, t in 1:mData.T], q2[i,t] >= pu2[i,t] * mData.pmin[i] + x[i,t] * (-mData.s[t] - mData.cc[i]) -
                                                                    (-mData.s[t] - mData.cc[i]) * mData.pmin[i] * u[clusterRev[i],t]);
    @constraint(mp, q2Cons2[i in 1:n, t in 1:mData.T], q2[i,t] >= pu2[i,t] * mData.pmax[i]);
    @constraint(mp, q2Cons3[i in 1:n, t in 1:mData.T], q2[i,t] <= pu2[i,t] * mData.pmin[i]);
    @constraint(mp, q2Cons4[i in 1:n, t in 1:mData.T], q2[i,t] <= pu2[i,t] * mData.pmax[i] + x[i,t] * (-mData.s[t] - mData.cc[i]) -
                                                                    (-mData.s[t] - mData.cc[i]) * mData.pmax[i] * u[clusterRev[i],t]);

    @objective(mp, Min, sum(V[t] for t in 1:mData.T));
    return mp;
end

function mccormick_opt(mData, n, clusterList, clusterRev, Gamma, time_limit = 3600)
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

    # π*u variables
    @variable(mp, pu1[1:n, 1:mData.T] >= 0);
    @variable(mp, pu2[1:n, 1:mData.T] <= 0);

    # mccormick relaxation terms
    @variable(mp, q1[1:n, 1:mData.T]);
    @variable(mp, q2[1:n, 1:mData.T]);

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

    # bilinear reformulation
    @constraint(mp, zplusCons[t in 1:mData.T], zplus[t] >= sum(q1[i,t] * mData.beta[i] for i in 1:n) +
                        π4[t,1] * Gamma[t] * sum(mData.pmax[i] for i in 1:n));
    @constraint(mp, zminusCons[t in 1:mData.T], zminus[t] >= sum(q2[i,t] * mData.alpha[i] for i in 1:n) +
                        π4[t,2] * Gamma[t] * sum(mData.pmax[i] for i in 1:n));

    # dual feasibility constraints
    @constraint(mp, dualCons1_plus[i in 1:n, t in 1:mData.T], mData.h[t] - mData.cc[i] - π1[i,t] <= π4[t,1]);
    @constraint(mp, dualCons2_plus[i in 1:n, t in 1:mData.T], -(mData.h[t] - mData.cc[i]) + π1[i,t] <= π4[t,1]);
    @constraint(mp, dualCons1_minus[i in 1:n, t in 1:mData.T], (-mData.s[t] - mData.cc[i]) - π2[i,t] <= π4[t,2]);
    @constraint(mp, dualCons2_minus[i in 1:n, t in 1:mData.T], mData.s[t] + mData.cc[i] + π2[i,t] <= π4[t,2]);

    # mccormick relaxation
    πbar1 = Dict();
    πubar2 = Dict();
    for i in 1:n
        for t in 1:mData.T
            πbar1[i,t] = mData.h[t] - mData.cc[i];
            πubar2[i,t] = -mData.s[t] - mData.cc[i];
        end
    end

    # π*u reformulation
    @constraint(mp, pu1Cons1[i in 1:n, t in 1:mData.T], pu1[i,t] <= π1[i,t]);
    @constraint(mp, pu1Cons2[i in 1:n, t in 1:mData.T], pu1[i,t] <= (mData.h[t] - mData.cc[i]) * u[clusterRev[i],t]);
    @constraint(mp, pu1Cons3[i in 1:n, t in 1:mData.T], pu1[i,t] >= π1[i,t] - (mData.h[t] - mData.cc[i]) * (1 - u[clusterRev[i],t]));

    @constraint(mp, pu2Cons1[i in 1:n, t in 1:mData.T], pu2[i,t] >= π2[i,t]);
    @constraint(mp, pu2Cons2[i in 1:n, t in 1:mData.T], pu2[i,t] >= (-mData.s[t] - mData.cc[i]) * u[clusterRev[i],t]);
    @constraint(mp, pu2Cons3[i in 1:n, t in 1:mData.T], pu2[i,t] <= π2[i,t] - (-mData.s[t] - mData.cc[i]) * (1 - u[clusterRev[i],t]));

    # McCormick relaxationn
    @constraint(mp, q1Cons1[i in 1:n, t in 1:mData.T], q1[i,t] >= pu1[i,t] * mData.pmin[i]);
    @constraint(mp, q1Cons2[i in 1:n, t in 1:mData.T], q1[i,t] >= pu1[i,t] * mData.pmax[i] + x[i,t] * (mData.h[t] - mData.cc[i]) -
                                                                    (mData.h[t] - mData.cc[i]) * mData.pmax[i] * u[clusterRev[i],t]);
    @constraint(mp, q1Cons3[i in 1:n, t in 1:mData.T], q1[i,t] <= pu1[i,t] * mData.pmin[i] + x[i,t] * (mData.h[t] - mData.cc[i]) -
                                                                    (mData.h[t] - mData.cc[i]) * mData.pmin[i] * u[clusterRev[i],t]);
    @constraint(mp, q1Cons4[i in 1:n, t in 1:mData.T], q1[i,t] <= pu1[i,t] * mData.pmax[i]);

    @constraint(mp, q2Cons1[i in 1:n, t in 1:mData.T], q2[i,t] >= pu2[i,t] * mData.pmin[i] + x[i,t] * (-mData.s[t] - mData.cc[i]) -
                                                                    (-mData.s[t] - mData.cc[i]) * mData.pmin[i] * u[clusterRev[i],t]);
    @constraint(mp, q2Cons2[i in 1:n, t in 1:mData.T], q2[i,t] >= pu2[i,t] * mData.pmax[i]);
    @constraint(mp, q2Cons3[i in 1:n, t in 1:mData.T], q2[i,t] <= pu2[i,t] * mData.pmin[i]);
    @constraint(mp, q2Cons4[i in 1:n, t in 1:mData.T], q2[i,t] <= pu2[i,t] * mData.pmax[i] + x[i,t] * (-mData.s[t] - mData.cc[i]) -
                                                                    (-mData.s[t] - mData.cc[i]) * mData.pmax[i] * u[clusterRev[i],t]);

    @objective(mp, Min, sum(V[t] for t in 1:mData.T));
    return mp;
end

function solve_mccormick(mData, n, clusterList, clusterRev, Gamma, relax_mode = 1, time_limit = 3600)
    # solve the McCormick relaxation
    if relax_mode == 1
        mcp = mccormick_opt_relax(mData, n, clusterList, clusterRev, Gamma, time_limit);
    elseif relax_mode == 2
        mcp = mccormick_opt(mData, n, clusterList, clusterRev, Gamma, time_limit);
    end

    start_time = time();
    optimize!(mcp);
    elapsed_time = time() - start_time;

    # return status, solution, and optimal value
    mcp_status = termination_status(mcp);
    mcp_lb = objective_bound(mcp);
    mcp_ub = objective_value(mcp);
    x_val = value.(mcp[:x]);
    π1_val = value.(mcp[:π1]);
    π2_val = value.(mcp[:π2]);
    π4_val = value.(mcp[:π4]);
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
    return mcp_status, elapsed_time, mcp_lb, mcp_ub, xDict, πDict;
end
