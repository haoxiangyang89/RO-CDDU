function [obj,Pi,Lam,Mu,Tht,ImprUppTime]=ImprDual(myData,N,T,p,Zeta)
try
    %Construct the solution
    Pi = zeros(N,T);
    Lam = zeros(N,T);
    Mu = zeros(T,1);
    Tht = zeros(T,1);
    %Sorted value  
    %Construct the parameter for later model
    SortedH = zeros(T,N+1);
    SortedS = zeros(T,N+1);
    gamH = zeros(T,N);
    etaS = zeros(T,N);
    %Sorted index
    phi = zeros(T,N);
    psi = zeros(T,N);
    for t = 1:T
        myH = abs(myData.('HoldingCost')(t) - myData.('CommitmentCost'));
        [SortedH(t,1:N),phi(t,:)] = sort(myH,'descend');
        myS = abs(myData.('ShortageCost')(t)+myData.('CommitmentCost'));
        [SortedS(t,1:N),psi(t,:)] = sort(myS,'descend');
        ck = (myData.('HoldingCost')(t)<= myData.('CommitmentCost'));
        gamH(t,:) = abs(myData.('alpha')).*ck + abs(myData.('beta')).*(1-ck);
        ck = (-myData.('ShortageCost')(t)<=myData.('CommitmentCost'));
        etaS(t,:) = abs(myData.('alpha')).*ck + abs(myData.('beta')).*(1-ck);
    end
    tic;
    Y = zeros(T,1);
    for t = 1:T
        tempPi = zeros(N,N+1);
        %Subproblem H 
        tempY1 = zeros(N+1,1);
        for k = 1:N+1
            tempu = SortedH(t,k);
            for i = 1:N
                tempPi(phi(t,i),k) = max(SortedH(t,i) - tempu,0);
                tempY1(k) = tempY1(k)+tempPi(phi(t,i),k)*gamH(t,phi(t,i))*p(phi(t,i),t);
            end
            tempY1(k) = tempY1(k)+ tempu*Zeta;
        end
        [Y1,minIdx] = min(tempY1);
        Pi(:,t) = tempPi(:,minIdx);
        Mu(t) = SortedH(t,minIdx);
        
        
        %Subproblem S
        tempY2 = zeros(N+1,1);
        tempPi = zeros(N,N+1);
        for k = 1:N+1
            tempu = SortedS(t,k);
            for i = 1:N
                tempPi(psi(t,i),k) = max(SortedS(t,i) - tempu,0);
                tempY2(k) = tempY2(k)+tempPi(psi(t,i),k)*etaS(t,psi(t,i))*p(psi(t,i),t);
            end
            tempY2(k) = tempY2(k)+ tempu*Zeta;
        end
        [Y2,minIdx] = min(tempY2);
        Lam(:,t) = tempPi(:,minIdx);
        Tht(t) = SortedS(t,minIdx);
        
        Y(t) = max(Y1+sum(myData.('HoldingCost')(t)*p(1:N,t))-...
            myData.('CommitmentCost')(1:N)*p(1:N,t)-...
            myData.('HoldingCost')(t)*myData.('D')(t),...
            Y2-sum(myData.('ShortageCost')(t)*p(1:N,t))...
            -myData.('CommitmentCost')(1:N)*p(1:N,t)+...
            myData.('ShortageCost')(t)*myData.('D')(t));
    end
    ImprUppTime = toc;
    obj = sum(Y);
catch err
    throw(err);
end
end