# original bilinear formulation
function bilinear_opt(mData, n, clusterList, clusterRev, Gamma, time_limit = 3600)
    mp = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(GUROBI_ENV),
                                         "Nonconvex" => 2,
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
    @variable(mp, π1[1:n, 1:mData.T, 1:2] >= 0);
    @variable(mp, π2[1:n, 1:mData.T, 1:2] <= 0);
    @variable(mp, π3[1:n, 1:mData.T, 1:2]);
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

    # bilinear reformulation
    @constraint(mp, zplusCons[t in 1:mData.T], zplus[t] >= sum(π1[i,t,1] * mData.beta[i] * x[i,t] + π2[i,t,1] * mData.alpha[i] * x[i,t] for i in 1:n) +
                        π4[t,1] * Gamma[t] * sum(mData.pmax[i] for i in 1:n));
    @constraint(mp, zminusCons[t in 1:mData.T], zminus[t] >= sum(π1[i,t,2] * mData.beta[i] * x[i,t] + π2[i,t,2] * mData.alpha[i] * x[i,t] for i in 1:n) +
                        π4[t,2] * Gamma[t] * sum(mData.pmax[i] for i in 1:n));

    # dual feasibility constraints
    @constraint(mp, dualCons1_plus[i in 1:n, t in 1:mData.T], -(π1[i,t,1] + π2[i,t,1] + π3[i,t,1]) + mData.h[t] - mData.cc[i] == 0);
    @constraint(mp, dualCons1_minus[i in 1:n, t in 1:mData.T], -(π1[i,t,2] + π2[i,t,2] + π3[i,t,2]) - mData.s[t] - mData.cc[i] == 0);
    @constraint(mp, dualCons2[i in 1:n, t in 1:mData.T, l in 1:2], π4[t,l] >= π3[i,t,l]);
    @constraint(mp, dualCons3[i in 1:n, t in 1:mData.T, l in 1:2], π4[t,l] >= -π3[i,t,l]);

    @objective(mp, Min, sum(V[t] for t in 1:mData.T));
    return mp;
end

function bilinear_opt_relax(mData, n, clusterList, clusterRev, Gamma, time_limit = 3600)
    mp = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(GUROBI_ENV),
                                         "Nonconvex" => 2,
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
    @variable(mp, π1[1:n, 1:mData.T, 1:2] >= 0);
    @variable(mp, π2[1:n, 1:mData.T, 1:2] <= 0);
    @variable(mp, π3[1:n, 1:mData.T, 1:2]);
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

    # bilinear reformulation
    @constraint(mp, zplusCons[t in 1:mData.T], zplus[t] >= sum(π1[i,t,1] * mData.beta[i] * x[i,t] + π2[i,t,1] * mData.alpha[i] * x[i,t] for i in 1:n) +
                        π4[t,1] * Gamma[t] * sum(mData.pmax[i] for i in 1:n));
    @constraint(mp, zminusCons[t in 1:mData.T], zminus[t] >= sum(π1[i,t,2] * mData.beta[i] * x[i,t] + π2[i,t,2] * mData.alpha[i] * x[i,t] for i in 1:n) +
                        π4[t,2] * Gamma[t] * sum(mData.pmax[i] for i in 1:n));

    # dual feasibility constraints
    @constraint(mp, dualCons1_plus[i in 1:n, t in 1:mData.T], -(π1[i,t,1] + π2[i,t,1] + π3[i,t,1]) + mData.h[t] - mData.cc[i] == 0);
    @constraint(mp, dualCons1_minus[i in 1:n, t in 1:mData.T], -(π1[i,t,2] + π2[i,t,2] + π3[i,t,2]) - mData.s[t] - mData.cc[i] == 0);
    @constraint(mp, dualCons2[i in 1:n, t in 1:mData.T, l in 1:2], π4[t,l] <= π3[i,t,l]);
    @constraint(mp, dualCons3[i in 1:n, t in 1:mData.T, l in 1:2], π4[t,l] <= -π3[i,t,l]);

    @objective(mp, Min, sum(V[t] for t in 1:mData.T));
    return mp;
end

function bilinear_simple(mData, n, clusterList, clusterRev, Gamma, time_limit = 3600)
    mp = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(GUROBI_ENV),
                                         "Nonconvex" => 2,
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

    # bilinear reformulation
    @constraint(mp, zplusCons[t in 1:mData.T], zplus[t] >= sum(π1[i,t] * mData.beta[i] * x[i,t] for i in 1:n) +
                        π4[t,1] * Gamma[t] * sum(mData.pmax[i] for i in 1:n));
    @constraint(mp, zminusCons[t in 1:mData.T], zminus[t] >= sum(π2[i,t] * mData.alpha[i] * x[i,t] for i in 1:n) +
                        π4[t,2] * Gamma[t] * sum(mData.pmax[i] for i in 1:n));

    # dual feasibility constraints
    @constraint(mp, dualCons1_plus[i in 1:n, t in 1:mData.T], mData.h[t] - mData.cc[i] - π1[i,t] <= π4[t,1]);
    @constraint(mp, dualCons2_plus[i in 1:n, t in 1:mData.T], -(mData.h[t] - mData.cc[i]) + π1[i,t] <= π4[t,1]);
    @constraint(mp, dualCons1_minus[i in 1:n, t in 1:mData.T], (- mData.s[t] - mData.cc[i]) - π2[i,t] <= π4[t,2]);
    @constraint(mp, dualCons2_minus[i in 1:n, t in 1:mData.T], mData.s[t] + mData.cc[i] + π2[i,t] <= π4[t,2]);

    @objective(mp, Min, sum(V[t] for t in 1:mData.T));
    return mp;
end

function solve_bilinear(mData, n, clusterList, clusterRev, Gamma, time_limit)
    # create/solve a bilinear program
    bps = bilinear_simple(mData, n, clusterList, clusterRev, Gamma, time_limit);
    start_time = time();
    optimize!(bps);
    elapsed_time = time() - start_time;

    # return status, solution, and optimal value
    bps_status = termination_status(bps);
    bps_lb = objective_bound(bps);
    bps_ub = objective_value(bps);
    x_val = value.(bps[:x]);
    π1_val = value.(bps[:π1]);
    π2_val = value.(bps[:π2]);
    π4_val = value.(bps[:π4]);
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
    return bps_status, elapsed_time, bps_lb, bps_ub, xDict, πDict;
end
