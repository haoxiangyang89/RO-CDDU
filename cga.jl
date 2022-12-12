# column generation
function cga_ini_xprob(mData, n, clusterList, clusterRev, Gamma)
    # create the initial x problem for column generation
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

    @objective(mp, Min, sum(V[t] for t in 1:mData.T));

    return mp;
end

function set_ini(mData, xprob_set, n, xIni, μIni)
    for t in 1:mData.T
        for i in 1:n
            set_start_value(xprob_set[:x][i,t], xIni[i,t]);
        end
    end
    for item in μIni
        t = item[1];
        cs = item[2];
        l = item[3];
        set_start_value(xprob_set[:μ][t,cs,l], 1);
    end
    return xprob_set;
end

function solve_x(mData, xprob, πList, n, clusterList, clusterRev, Gamma, M1, M2, return_prob = 0, time_limit = 10800, xIni = Dict(), μIni = [])
    xprob_temp = copy(xprob);
    set_optimizer(xprob_temp, () -> Gurobi.Optimizer(GUROBI_ENV));
    set_optimizer_attribute(xprob_temp, "TimeLimit", time_limit);
    set_optimizer_attribute(xprob_temp, "Threads", 30);

    s = length(πList);

    # incomplete column constraints
    # bilinear reformulation
    @variable(xprob_temp, μ[1:mData.T, 1:s, 1:2], Bin);
    @constraint(xprob_temp, μCons[t in 1:mData.T, l in 1:2], sum(μ[t, cs, l] for cs in 1:s) == 1);

    # add the current columns
    @constraint(xprob_temp, zplusCons[t in 1:mData.T, cs in 1:s], xprob_temp[:zplus][t] >= sum(πList[cs][1][i,t] * mData.beta[i] * xprob_temp[:x][i,t] for i in 1:n) +
                        πList[cs][4][t,1] * Gamma[t] * sum(mData.pmax[i] for i in 1:n) - M1[t] * (1 - μ[t,cs,1]));
    @constraint(xprob_temp, zminusCons[t in 1:mData.T, cs in 1:s], xprob_temp[:zminus][t] >= sum(πList[cs][2][i,t] * mData.alpha[i] * xprob_temp[:x][i,t] for i in 1:n) +
                        πList[cs][4][t,2] * Gamma[t] * sum(mData.pmax[i] for i in 1:n) - M2[t] * (1 - μ[t,cs,2]));

    if return_prob == 0
        if xIni != Dict()
            xprob_temp = set_ini(mData, xprob_temp, n, xIni, μIni);
        end
        optimize!(xprob_temp);
        xhat = Dict();
        for i in 1:n
            for t in 1:mData.T
                xhat[i,t] = value(xprob_temp[:x][i,t]);
            end
        end

        μhat = [];
        μvalue = value.(xprob_temp[:μ]);
        for t in 1:mData.T
            for cs in 1:s
                for l in 1:2
                    if abs(μvalue[t,cs,l] - 1) <= 1e-4
                        push!(μhat, (t,cs,l));
                    end
                end
            end
        end

        return xhat, objective_value(xprob_temp), μhat;
    else
        return xprob_temp;
    end
end

function solve_x_new(mData, xprob, πList, n, clusterList, clusterRev, Gamma, return_prob = 0, time_limit = 10800)
    xprob_temp = copy(xprob);
    set_optimizer(xprob_temp, () -> Gurobi.Optimizer(GUROBI_ENV));
    set_optimizer_attribute(xprob_temp, "TimeLimit", time_limit);
    set_optimizer_attribute(xprob_temp, "Threads", 30);

    s = length(πList);

    # incomplete column constraints
    # bilinear reformulation
    @variable(xprob_temp, μ[1:mData.T, 1:s, 1:2], Bin);
    @variable(xprob_temp, γ[i in 1:n, 1:mData.T, 1:s, 1:2] >= 0);
    @constraint(xprob_temp, μCons[t in 1:mData.T, l in 1:2], sum(μ[t, cs, l] for cs in 1:s) == 1);

    # add the current columns
    @constraint(xprob_temp, zplusCons[t in 1:mData.T], xprob_temp[:zplus][t] >= sum(sum(πList[cs][1][i,t] * mData.beta[i] * xprob_temp[:γ][i, t, cs, 1] for i in 1:n) +
                        πList[cs][4][t,1] * μ[t, cs, 1] * Gamma[t] * sum(mData.pmax[i] for i in 1:n) for cs in 1:s));
    @constraint(xprob_temp, zminusCons[t in 1:mData.T], xprob_temp[:zminus][t] >= sum(sum(πList[cs][2][i,t] * mData.alpha[i] * xprob_temp[:γ][i, t, cs, 2] for i in 1:n) +
                        πList[cs][4][t,2] * μ[t, cs, 2] * Gamma[t] * sum(mData.pmax[i] for i in 1:n) for cs in 1:s));

    # add the linearization constraints
    @constraint(xprob_temp, gammaCon1[i in 1:n, t in 1:mData.T, cs in 1:s, l in 1:2], xprob_temp[:γ][i, t, cs, l] <= xprob_temp[:μ][t, cs, l] * mData.pmax[i]);
    @constraint(xprob_temp, gammaCon2[i in 1:n, t in 1:mData.T, cs in 1:s, l in 1:2], xprob_temp[:γ][i, t, cs, l] <= xprob_temp[:x][i, t]);
    @constraint(xprob_temp, gammaCon3[i in 1:n, t in 1:mData.T, cs in 1:s, l in 1:2], xprob_temp[:γ][i, t, cs, l] >= xprob_temp[:x][i, t] - (1 - xprob_temp[:μ][t, cs, l]) * mData.pmax[i]);

    if return_prob == 0
        optimize!(xprob_temp);
        xhat = Dict();
        for i in 1:n
            for t in 1:mData.T
                xhat[i,t] = value(xprob_temp[:x][i,t]);
            end
        end

        return xhat, objective_value(xprob_temp);
    else
        return xprob_temp;
    end
end


function solve_pi(mData, xhat, n, clusterList, clusterRev, Gamma)
    # given an x, solve for the pi
    xmax = sum(mData.pmax[i] for i in 1:n);
    π1hat = Dict();
    π2hat = Dict();
    π4hat = Dict();
    zplushat = Dict();
    zminushat = Dict();

    for t in 1:mData.T
        hcList = sort(unique([mData.h[t] - mData.cc[i] for i in 1:n]));
        pushfirst!(hcList,0.0);
        stopBool = false;
        zvalue = Inf;
        hcInd = 1;
        πbest = [[0.0],0.0];
        while (!(stopBool))&(hcInd <= length(hcList))
            π4 = hcList[hcInd];
            π1 = [max(mData.h[t] - mData.cc[i] - π4, 0) for i in 1:n];
            zvalue_new = sum(π1[i] * mData.beta[i] * xhat[i,t] for i in 1:n) + π4 * Gamma[t] * xmax;
            if zvalue_new < zvalue
                # update the π value
                zvalue = zvalue_new;
                hcInd += 1;
                πbest[1] = π1;
                πbest[2] = π4;
            else
                stopBool = true;
                hcInd -= 1;
            end
        end
        π4hat[t,1] = πbest[2];
        for i in 1:n
            π1hat[i,t] = πbest[1][i];
        end
        zplushat[t] = zvalue;
    end

    for t in 1:mData.T
        scList = sort(unique([mData.s[t] + mData.cc[i] for i in 1:n]));
        pushfirst!(scList,0.0);
        stopBool = false;
        zvalue = Inf;
        scInd = 1;
        πbest = [[0.0],0.0];
        while (!(stopBool))&(scInd <= length(scList))
            π4 = scList[scInd];
            π2 = [min(-mData.s[t] - mData.cc[i] + π4, 0) for i in 1:n];
            zvalue_new = sum(π2[i] * mData.alpha[i] * xhat[i,t] for i in 1:n) + π4 * Gamma[t] * xmax;
            if zvalue_new < zvalue
                # update the π value
                zvalue = zvalue_new;
                scInd += 1;
                πbest[1] = π2;
                πbest[2] = π4;
            else
                stopBool = true;
                scInd -= 1;
            end
        end
        π4hat[t,2] = πbest[2];
        for i in 1:n
            π2hat[i,t] = πbest[1][i];
        end
        zminushat[t] = zvalue;
    end

    Vhat = Dict();
    for t in 1:mData.T
        Vhat[t] = max(sum((mData.h[t] - mData.cc[i]) * xhat[i,t] for i in 1:n) - mData.h[t] * mData.D[t] + zplushat[t],
        sum((-mData.s[t] - mData.cc[i]) * xhat[i,t] for i in 1:n) + mData.s[t] * mData.D[t] + zminushat[t]);
    end

    πhat = Dict();
    πhat[1] = π1hat;
    πhat[2] = π2hat;
    πhat[4] = π4hat;

    return πhat;
end

function main_cga(mData, n, clusterList, clusterRev, Gamma, mode = 1, time_limit = 36000)
    # the process to run column generation

    # initialize a π
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
    if mode == 1
        xhat, vbest, μhat = solve_x(mData, xprob, πList, n, clusterList, clusterRev, Gamma, M1, M2);
    else
        xhat, vbest = solve_x_new(mData, xprob, πList, n, clusterList, clusterRev, Gamma);
    end
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

        if sum(sum(abs(xhat[i,t] - xhat_new[i,t]) for t in 1:mData.T) for i in 1:n) <= 1e-4
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

    return elapsed_time, vbest, xbest, πbest;
end
