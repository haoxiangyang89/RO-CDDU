# auxiliary functions

function genClusters(n, m)
    # split the n resource into m clusters uniformly
    counter = 0;
    cluster_counter = 1;
    clusterList = [i for i in 1:m];
    cluster_temp = [];
    cluster_size = Int64(round(n/m));
    clusterRev = Dict();
    clusterInfo = Dict();
    for i in 1:n
        counter += 1;
        push!(cluster_temp,i);
        clusterRev[i] = cluster_counter;
        if counter == cluster_size
            counter = 0;
            clusterInfo[cluster_counter] = cluster_temp;
            cluster_counter += 1;
            cluster_temp = [];
        end
    end
    return clusterList, clusterInfo, clusterRev;
end

function combineGen(mData)
    # combine the generator with the same cost
    uc = unique(mData.cc);

    for i in 1:length(uc)
        sameInd = findall(x -> x==uc[i], mData.cc);
        for j in 1:length(sameInd)
            mData.cc[sameInd[j]] = round(mData.cc[sameInd[j]] + 0.001*j, digits = 3);
        end
    end

    return mData;
end

function getM(n, mData, Gamma)
    M1 = Dict();
    M2 = Dict();
    for t in 1:mData.T
        hcList = sort(unique([mData.h[t] - mData.cc[i] for i in 1:n]));
        pushfirst!(hcList,0.0);
        scList = sort(unique([mData.s[t] + mData.cc[i] for i in 1:n]));
        pushfirst!(scList,0.0);
        M1best = -Inf;
        M2best = -Inf;

        for ihc in 1:length(hcList)
            π4 = hcList[ihc];
            π1 = [max(mData.h[t] - mData.cc[i] - π4, 0) for i in 1:n];
            zvalue_new = sum(π1[i] * mData.beta[i] * mData.pmax[i] + π4 * Gamma[t] * mData.pmax[i] for i in 1:n);
            if zvalue_new > M1best
                M1best = zvalue_new;
            end
        end
        for isc in 1:length(scList)
            π4 = scList[isc];
            π2 = [min(-mData.s[t] - mData.cc[i] + π4, 0) for i in 1:n];
            zvalue_new = sum(π2[i] * mData.alpha[i] * mData.pmax[i] + π4 * Gamma[t] * mData.pmax[i] for i in 1:n);
            if zvalue_new > M2best
                M2best = zvalue_new;
            end
        end

        M1[t] = M1best;
        M2[t] = M2best;
    end
    return M1, M2;
end

function scenGen(mData, n, m)
    # generate m scenarios of uniformly distributed loads
    scenDict = Dict();
    for j in 1:m
        # for each scenario, generate random var
        scenDict[j] = rand(n, mData.T);
    end
    return scenDict;
end
