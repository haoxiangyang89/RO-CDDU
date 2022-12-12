# deterministic optimization
# clusterList: list of clusters
# clusterInfo is a dictionary: cluster_no: [set of resource within the cluster]
# clusterRev is a dictionary: reverse of clusterInfo

function det_opt(mData, n, clusterList, clusterRev)
    mp = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(GUROBI_ENV)));
    # set up the variables:
    #   x: amount of demand response
    #   u: commitment indicator
    #   w,v: increase/decrease indicator
    @variable(mp, x[i in 1:n, t in 1:mData.T]);
    @variable(mp, zplus[t in 1:mData.T] >= 0);
    @variable(mp, zminus[t in 1:mData.T] >= 0);
    @variable(mp, u[i in clusterList, t in 1:mData.T], Bin);
    @variable(mp, w[i in clusterList, t in 2:mData.T], Bin);

    # commitment bounds
    @constraint(mp, commitment_up[i in 1:n, t in 1:mData.T], x[i,t] <= mData.pmax[i] * u[clusterRev[i],t]);
    @constraint(mp, commitment_down[i in 1:n, t in 1:mData.T], x[i,t] >= mData.pmin[i] * u[clusterRev[i],t]);

    # ramping constraints
    @constraint(mp, ramp_up[i in 1:n, t in 1:(mData.T - 1)], x[i,t + 1] - x[i,t] <= mData.RU[i] * u[clusterRev[i],t + 1]);
    @constraint(mp, ramp_down[i in 1:n, t in 1:(mData.T - 1)], x[i,t + 1] - x[i,t] >= -mData.RD[i] * u[clusterRev[i],t]);

    # cooling down constraints
    @constraint(mp, w1[i in 1:n, t in 2:mData.T], x[i,t] - x[i,t - 1] <= mData.RU[i] * w[clusterRev[i],t]);
    @constraint(mp, w2[i in 1:n, t in 2:mData.T], x[i,t] - x[i,t - 1] >= -mData.RD[i] * (1 - w[clusterRev[i],t]));
    @constraint(mp, w3[i in 1:n, t in 2:mData.T, τ in t:min(t + Int64(mData.upmin[i]) - 1,mData.T)], x[i,τ] - x[i,τ - 1] >= -mData.RD[i] * (1 - w[clusterRev[i],t]));
    @constraint(mp, w4[i in 1:n, t in 2:mData.T, τ in t:min(t + Int64(mData.downmin[i]) - 1,mData.T)], x[i,τ] - x[i,τ - 1] <= mData.RU[i] * w[clusterRev[i],t]);

    # tightening constraint
    # @constraint(mp, tighten1[i in 1:n, t in 2:mData.T], w[i,t] <= u[i,t]);
    # @constraint(mp, tighten2[i in 1:n, t in 2:mData.T], 1 - w[i,t] <= u[i,t - 1]);

    # absolute value of sum(x) - D
    @constraint(mp, costPos1[t in 1:mData.T], zplus[t] >= sum(x[i,t] for i in 1:n) - mData.D[t]);
    @constraint(mp, costPos2[t in 1:mData.T], zminus[t] >= -sum(x[i,t] for i in 1:n) + mData.D[t]);

    @objective(mp, Min, sum(mData.h[t] * zplus[t] + mData.s[t] * zminus[t] - sum(mData.cc[i] * x[i,t] for i in 1:n) for t in 1:mData.T));
    return mp;
end

function det_opt_2(mData, n, clusterList, clusterRev)
    mp = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(GUROBI_ENV)));
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

    # absolute value of sum(x) - D
    @constraint(mp, costPos1[t in 1:mData.T], zplus[t] >= sum(x[i,t] for i in 1:n) - mData.D[t]);
    @constraint(mp, costPos2[t in 1:mData.T], zminus[t] >= -sum(x[i,t] for i in 1:n) + mData.D[t]);

    @objective(mp, Min, sum(mData.h[t] * zplus[t] + mData.s[t] * zminus[t] - sum(mData.cc[i] * x[i,t] for i in 1:n) for t in 1:mData.T));
    return mp;
end

function eval_opt(mData, n, xHat, scenData)
    # given a solution and a random ξ, calculate the objective function
    xi = Dict();
    total_cost = 0;
    for t in 1:mData.T
        for i in 1:n
            xi[i,t] = scenData[i,t]*(mData.beta[i] - mData.alpha[i]) + mData.alpha[i];
            total_cost -= mData.cc[i] * (xHat[i,t]*(1 + xi[i,t]));
        end
        mismatch = sum(xHat[i,t]*(1 + xi[i,t]) for i in 1:n) - mData.D[t];
        if mismatch > 0
            total_cost += mData.h[t] * mismatch;
        else
            total_cost += mData.s[t] * (-mismatch);
        end
    end
    return total_cost;
end
