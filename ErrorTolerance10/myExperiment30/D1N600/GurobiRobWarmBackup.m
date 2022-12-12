function [obj,p,myRunTime]=GurobiRob(myData,N,T,Zeta,solvec)
    model.modelName = 'myDRModel';
    model.modelsense = 'Min';
    
    M=max(myData.('pmax'));
    %Construct the parameter for later model
    %Sorted value  
    SortedH = zeros(T,N);
    SortedS = zeros(T,N);
    gamH = zeros(T,N);
    etaS = zeros(T,N);
    %Sorted index
    phi = zeros(T,N);
    psi = zeros(T,N);
    for t = 1:T
        myH = abs(myData.('HoldingCost')(t) - myData.('CommitmentCost'));
        [SortedH(t,:),phi(t,:)] = sort(myH,'descend');
        myS = abs(myData.('ShortageCost')(t)+myData.('CommitmentCost'));
        [SortedS(t,:),psi(t,:)] = sort(myS,'descend');
        ck = (myData.('HoldingCost')(t)<= myData.('CommitmentCost'));
        gamH(t,:) = abs(myData.('alpha')).*ck + abs(myData.('beta')).*(1-ck);
        ck = (-myData.('ShortageCost')(t)>=myData.('CommitmentCost'));
        etaS(t,:) = abs(myData.('alpha')).*ck + abs(myData.('beta')).*(1-ck);
    end
     Mt=sum((max(myData.('ShortageCost'))+...
        max(myData.('CommitmentCost'))).*myData.('pmax'))*max(max(gamH));
    
    %Sequence of decision varialbe in matrix 
    %Y(t),p(i,t), w(i,t), v(i,t),  x(i,t), W(i,t), U(i,t), nu(i,t),mu(i,t)
    %prepare feasible solutions for W,U,nu,mu
    WUNM = zeros(4*(N+1)*T,1);
     for t=  1:T
         tempW = zeros(N+1,1);
         tempU = zeros(N+1,1);
         for k = 1:N
            for i = 1:k-1
                tempW(k) = tempW(k)+(SortedH(t,i) - SortedH(t,k))...
                    *gamH(t,phi(t,i))*solvec(T+(t-1)*N+phi(t,i));
                tempU(k) = tempU(k)+(SortedS(t,i) - SortedS(t,k))...
                    *etaS(t,psi(t,i))*solvec(T+(t-1)*N+psi(t,i));
            end
            tempW(k)=tempW(k)+SortedH(t,k)*Zeta;
            tempU(k)=tempU(k)+SortedS(t,k)*Zeta;
        end

        for i = 1:N
            tempW(N+1) = tempW(N+1) +...
                SortedH(t,i)*gamH(t,phi(t,i))*solvec(T+(t-1)*N+phi(t,i));
            tempU(N+1) = tempU(N+1)+...
                 SortedS(t,i)*etaS(t,psi(t,i))*solvec(T+(t-1)*N+psi(t,i));
                
        end
        [sortW,inxW] = sort(tempW,'ascend');
        [sortU,inxU] = sort(tempU,'ascend');
        %W(i,t)
        WUNM((t-1)*(N+1)+1:t*(N+1)) = max(tempW - Mt,0);
        %U(i,t)
        WUNM((T*(N+1)+(t-1)*(N+1)+1):(T*(N+1)+t*(N+1))) = max(tempU - Mt,0);
        WUNM((t-1)*(N+1)+inxW(1)) = sortW(1);
        WUNM(T*(N+1)+(t-1)*(N+1)+inxU(1)) = sortU(1);
        %nu(i,t),mu(i,t)
        WUNM(2*T*(N+1)+(t-1)*(N+1)+inxW(1)) = 1;
        WUNM(3*T*(N+1)+(t-1)*(N+1)+inxU(1)) = 1;
     end
     myStart = [solvec;WUNM];
    model.Start = myStart;
    
    myObjArr = [ones(T,1);zeros(4*N*T+4*(N+1)*T,1);]; 
    model.obj  = myObjArr;
    
    %Y(t),p(i,t), w(i,t), v(i,t),  x(i,t),W(i,t), U(i,t), nu(i,t),mu(i,t)
    %Set the lower bound for the decision variable
    DVLowerBound = [-inf*ones(T,1);zeros(4*N*T,1);zeros(2*(N+1)*T,1);zeros(2*(N+1)*T,1)];
    model.lb    = DVLowerBound;
    %Set the upper bound for the decision variables
    DVUpperBound = [inf*ones(T+N*T,1);ones(3*N*T,1);inf*ones(2*(N+1)*T,1);ones(2*(N+1)*T,1)];
    model.ub    = DVUpperBound;
    %Y(t),p(i,t), w(i,t), v(i,t),  x(i,t),W(i,t), U(i,t), nu(i,t),mu(i,t)
    %Set the datatype for the decision variables
    Datatype = [repmat('C',[1,T+N*T]),repmat('B',[1,3*N*T]),repmat('C',[1,2*(N+1)*T]),...
        repmat('B',[1,2*(N+1)*T])];
    model.vtype = Datatype;

    %Count the total number of constraints
    NumOfCol = T+4*N*T+4*(N+1)*T;
    NumOfRow=0;
    %Objective constrints, Capacity constraints and ramping constraints
    NumOfRow=NumOfRow+2*((N+1)*T+N*T+T+T);
    NumOfRow = NumOfRow+2*N*T+(T-1)*N +N*T;
    
    %Upupmin constraints
    NumOfRow = NumOfRow+N*T;
    for i=1:N
        for t=1:T
            for tau=t:min(t+myData.('Upupmin')(i)-1,T)
                NumOfRow=NumOfRow+1;
            end
        end
    end
    
    %Downdownmin constraints
    NumOfRow = NumOfRow+N*T;
    for i=1:N
        for t=1:T
            for tau=t:min(t+myData.('Downdownmin')(i)-1,T)
                NumOfRow=NumOfRow+1;
            end
        end
    end
    %Set the main constraint matrix.
    ConstraintBody=sparse(NumOfRow,NumOfCol);   

    
    %Row index
    ii = 0;
    %Y(t),p(i,t), w(i,t), v(i,t),  x(i,t),W(i,t), U(i,t), nu(i,t),mu(i,t)
    %Objective Constraint 1(1)
    for t = 1:T
        for k = 1:N+1
            ii=ii+1;
            ConstraintBody(ii,t) = 1;
            ConstraintBody(ii,T+4*N*T+(t-1)*(N+1)+k) = -1;
            ConstraintBody(ii,T+(t-1)*N+(1:N)) = -myData.('HoldingCost')(t)+myData.('CommitmentCost')(1:N);
%             for i = 1:N 
%                 ConstraintBody(ii,T+(t-1)*N+i) = -myData.('HoldingCost')(t)+myData.('CommitmentCost')(i);
%             end
        end
    end
    %Objective Constraint 1(2)
    %Y(t),p(i,t), w(i,t), v(i,t),  x(i,t),W(i,t), U(i,t), nu(i,t),mu(i,t)
    for t=  1:T
        for k = 1:N
            ii=ii+1;
            ConstraintBody(ii,T+4*N*T+(t-1)*(N+1)+k) = 1;
            ConstraintBody(ii,T+4*N*T+2*(N+1)*T+(t-1)*(N+1)+k) = -Mt;
            ConstraintBody(ii,T+(t-1)*N+(phi(t,1:k-1))) = (-SortedH(t,1:k-1) + SortedH(t,k))...
                    .*gamH(t,phi(t,1:k-1));
%             for i = 1:k-1
%                 ConstraintBody(ii,T+(t-1)*N+phi(t,i)) = (-SortedH(t,i) + SortedH(t,k))...
%                     *gamH(t,phi(t,i));
%             end
        end
        ii = ii+1;
        ConstraintBody(ii,T+4*N*T+t*(N+1)) = 1; 
        ConstraintBody(ii,T+4*N*T+2*(N+1)*T+t*(N+1)) = -Mt;
        ConstraintBody(ii,T+(t-1)*N+phi(t,1:N)) =  -SortedH(t,1:N).*gamH(t,phi(t,1:N));
%         for i = 1:N
%             ConstraintBody(ii,T+(t-1)*N+phi(t,i)) =  -SortedH(t,i)*gamH(t,phi(t,i));
%         end
    end
    %Objective Constraint 1(3)
    for t = 1:T
        ii=ii+1;
        ConstraintBody(ii,T+4*N*T+2*(N+1)*T+(t-1)*(N+1)+(1:N+1)) = 1;
%         for i = 1:N+1
%             ConstraintBody(ii,T+4*N*T+2*(N+1)*T+(t-1)*(N+1)+i) = 1; 
%         end
    end
    
    %Y(t),p(i,t), w(i,t), v(i,t),  x(i,t),W(i,t), U(i,t), nu(i,t),mu(i,t)
    %Objective Constraint 2(1)
    for t = 1:T
        for k = 1:N+1
            ii=ii+1;
            ConstraintBody(ii,t) = 1;
            ConstraintBody(ii,T+4*N*T+(N+1)*T+(t-1)*(N+1)+k) = -1;
            ConstraintBody(ii,T+(t-1)*N+(1:N)) = myData.('ShortageCost')(t)+myData.('CommitmentCost')(1:N);
%             for i = 1:N 
%                 ConstraintBody(ii,T+(t-1)*N+i) = myData.('ShortageCost')(t)+myData.('CommitmentCost')(i);
%             end
        end
    end
    %Objective Constraint 2(2)
    %Y(t),p(i,t), w(i,t), v(i,t),  x(i,t),W(i,t), U(i,t), nu(i,t),mu(i,t)
    for t=  1:T
        for k = 1:N
            ii=ii+1;
            ConstraintBody(ii,T+4*N*T+(N+1)*T+(t-1)*(N+1)+k) = 1;
            ConstraintBody(ii,T+4*N*T+3*(N+1)*T+(t-1)*(N+1)+k) = -Mt;
            ConstraintBody(ii,T+(t-1)*N+psi(t,1:k-1)) = (-SortedS(t,1:k-1) + SortedS(t,k))...
                    .*etaS(t,psi(t,1:k-1));
%             for i = 1:k-1
%                 ConstraintBody(ii,T+(t-1)*N+psi(t,i)) = (-SortedS(t,i) + SortedS(t,k))...
%                     *etaS(t,psi(t,i));
%             end
        end
        ii = ii+1;
        ConstraintBody(ii,T+4*N*T+(N+1)*T+t*(N+1)) = 1; 
        ConstraintBody(ii,T+4*N*T+3*(N+1)*T+t*(N+1)) = -Mt;
        ConstraintBody(ii,T+(t-1)*N+psi(t,1:N)) = -SortedS(t,1:N).*etaS(t,psi(t,1:N));
%         for i = 1:N
%             ConstraintBody(ii,T+(t-1)*N+psi(t,i)) = -SortedS(t,i)*etaS(t,psi(t,i));
%         end
    end

    %Objective Constraint 2(3)
    %Y(t),p(i,t), w(i,t), v(i,t),  x(i,t),W(i,t), U(i,t), nu(i,t),mu(i,t)
    for t = 1:T
        ii=ii+1;
        ConstraintBody(ii,T+4*N*T+3*(N+1)*T+(t-1)*(N+1)+(1:N+1)) = 1;
%         for i = 1:N+1
%             ConstraintBody(ii,T+4*N*T+3*(N+1)*T+(t-1)*(N+1)+i) = 1;
%         end
    end
     %Y(t), p(i,t), w(i,t), v(i,t),  x(i,t),W(i,t), U(i,t), nu(i,t),mu(i,t)
    %pmin
    for i =1:N
        for t=1:T
            ii=ii+1;
            ConstraintBody(ii,T+(t-1)*N+i)=1;
            ConstraintBody(ii,T+3*N*T+(t-1)*N+i)=-myData.('pmin')(i);
        end
    end

     %Y(t), p(i,t), w(i,t), v(i,t),  x(i,t),W(i,t), U(i,t), nu(i,t),mu(i,t)
    %pmax
    for i =1:N
        for t=1:T
            ii=ii+1;
            ConstraintBody(ii,T+(t-1)*N+i)=-1;
            ConstraintBody(ii,T+3*N*T+(t-1)*N+i)=myData.('pmax')(i);
        end
    end
    %Y(t), p(i,t), w(i,t), v(i,t),  x(i,t),W(i,t), U(i,t), nu(i,t),mu(i,t)
    %RD
    for i = 1:N
        for t= 1:T-1
            ii=ii+1;
            ConstraintBody(ii,T+t*N+i)=1;
            ConstraintBody(ii,T+(t-1)*N+i)=-1;
            ConstraintBody(ii,T+3*N*T+(t-1)*N+i) = myData.('RD')(i);
        end
    end
    %Y(t), p(i,t), w(i,t), v(i,t),  x(i,t),W(i,t), U(i,t), nu(i,t),mu(i,t)
    %RU
    for i = 1:N
        ii=ii+1;
        ConstraintBody(ii,T+i)=-1;
        ConstraintBody(ii,T+3*N*T+i)=myData.('RU')(i);
        for t=2:T
            ii=ii+1;
            ConstraintBody(ii,T+(t-1)*N+i)=-1;
            ConstraintBody(ii,T+(t-2)*N+i)=1;
            ConstraintBody(ii,T+3*N*T+(t-1)*N+i)=myData.('RU')(i);
        end
    end
    
    %Y(t), p(i,t), w(i,t), v(i,t),  x(i,t),W(i,t), U(i,t), nu(i,t),mu(i,t)
    %Upupmin
    for i= 1:N
        ii=ii+1;
        ConstraintBody(ii,T+i) = -1;
        ConstraintBody(ii,T+N*T+i) = M;
        for t=2:T
            ii=ii+1;
            ConstraintBody(ii,T+(t-1)*N+i)=-1;
            ConstraintBody(ii,T+(t-2)*N+i)=1;
            ConstraintBody(ii,T+N*T+(t-1)*N+i)=M;
        end
    end
    for i=1:N
        ii=ii+1;
        ConstraintBody(ii,T+i)=1;
        ConstraintBody(ii,T+N*T+i)=-M;
        for tau = 2:min(myData.('Upupmin')(i),T)
            ii=ii+1;
            ConstraintBody(ii,T+(tau-1)*N+i)=1;
            ConstraintBody(ii,T+(tau-2)*N+i)=-1;
            ConstraintBody(ii,T+N*T+i)=-M;
        end
        for t= 2: T 
            for tau = t:min(t+myData.('Upupmin')(i)-1,T)
                ii=ii+1;
                ConstraintBody(ii,T+(tau-1)*N+i)=1;
                ConstraintBody(ii,T+(tau-2)*N+i)=-1;
                ConstraintBody(ii,T+N*T+(t-1)*N+i)=-M;
            end
        end
    end

    %Y(t), p(i,t), w(i,t), v(i,t),  x(i,t),W(i,t), U(i,t), nu(i,t),mu(i,t)
    %Downdownmin
    for i= 1:N
        ii=ii+1;
        ConstraintBody(ii,T+i) = 1;
        ConstraintBody(ii,T+2*N*T+i) = M;
        for t=2:T
            ii=ii+1;
            ConstraintBody(ii,T+(t-1)*N+i)=1;
            ConstraintBody(ii,T+(t-2)*N+i)=-1;
            ConstraintBody(ii,T+2*N*T+(t-1)*N+i)=M;
        end
    end
    for i =1:N
        ii=ii+1;
        ConstraintBody(ii,T+i)=-1;
        ConstraintBody(ii,T+2*N*T+i)=-M;
        for tau = 2:min(myData.('Downdownmin')(i),T)
            ii=ii+1;
            ConstraintBody(ii,T+(tau-1)*N+i)=-1;
            ConstraintBody(ii,T+(tau-2)*N+i)=1;
            ConstraintBody(ii,T+2*N*T+i)= -M;
        end
        for t = 2:T
            for tau = t:min(t+myData.('Downdownmin')(i)-1,T)
                ii=ii+1;
                ConstraintBody(ii,T+(tau-1)*N+i)=-1;
                ConstraintBody(ii,T+(tau-2)*N+i)=1;
                ConstraintBody(ii,T+2*N*T+(t-1)*N+i)=-M;
            end
        end
    end
    model.A =  ConstraintBody;
    
    
    %Constraintbody's bounds on A
    Alower=zeros(NumOfRow,1);
    
    ii=0;
    %Objective Constraint 1(1) 
    %Y(t),p(i,t), w(i,t), v(i,t),  x(i,t),W(i,t), U(i,t), nu(i,t),mu(i,t)
    for t = 1:T
        for k = 1:N+1
            ii=ii+1;
            Alower(ii,1)=-myData.('HoldingCost')(t)*myData.('D')(t); 
        end
    end
    %Objective Constraint 1(2)
    %Y(t),p(i,t), w(i,t), v(i,t),  x(i,t),W(i,t), U(i,t), nu(i,t),mu(i,t)
    for t=  1:T
        for k = 1:N
            ii=ii+1;
            Alower(ii,1) = SortedH(t,k)*Zeta - Mt;
        end
        ii = ii+1;
        Alower(ii,1) = -Mt;
    end
    %Objective Constraint 1(3)
    for t = 1:T
        ii=ii+1;
        Alower(ii,1) = 1;
    end
    
    %Objective Constraint 2(1)
    for t = 1:T
        for k = 1:N+1
            ii=ii+1;
            Alower(ii,1)=myData.('ShortageCost')(t)*myData.('D')(t); 
        end
    end
    %Objective Constraint 2(2)
    %Y(t),p(i,t), w(i,t), v(i,t),  x(i,t),W(i,t), U(i,t), nu(i,t),mu(i,t)
    for t=  1:T
        for k = 1:N
            ii=ii+1;
            Alower(ii,1) = SortedS(t,k)*Zeta - Mt;
        end
        ii = ii+1;
        Alower(ii,1) = -Mt;
    end
    %Objective Constraint 1(3)
    for t = 1:T
        ii=ii+1;
        Alower(ii,1) = 1;
    end
    
    
    
    %pmin
    for i =1:N
        for t=1:T
            ii=ii+1;
            Alower(ii,1)=0;
        end
    end
    %pmax
    for i =1:N
        for t=1:T
            ii=ii+1;
            Alower(ii,1)=0;
        end
    end
    %RD
    for i =1:N
        for t= 1:T-1
            ii=ii+1;
            Alower(ii,1)=0;
        end
    end
    %RU
    for i=1:N
        for t= 1:T
            ii=ii+1;
            Alower(ii,1)=0;
        end
    end
    
    %Upupmin
    for i = 1:N
        for t= 1:T
            ii=ii+1;
            Alower(ii,1)=0;
        end
    end
    for i =1:N
        for t= 1:T
            for tau= t: min(t+myData.('Upupmin')(i)-1,T)
                ii=ii+1;
                Alower(ii,1)=-M;
            end
        end
    end
    %Downdownmin
    for i =1:N
        for t=1:T
            ii=ii+1;
            Alower(ii,1)=0;
        end
    end
    for i=1:N
        for t=1:T
            for tau=t:min(t+myData.('Downdownmin')(i)-1,T)
                ii=ii+1;
                Alower(ii,1)=-M;
            end
        end
    end
    
    
    %**********************************
    %Testing feasiblity*****************
%     model.lb =myStart-10^(-5);
%     model.ub = myStart+10^(-5);
    %************************************

    model.rhs = Alower;
    model.sense = repmat('>',NumOfRow,1);
    
    clear params;
    params.MIPGap = 0.1;

    
    
    tic
    results = gurobi(model,params);
    myRunTime = toc;
    
    %Y(t), p(i,t), w(i,t), v(i,t),  x(i,t),W(i,t), U(i,t), nu(i,t),mu(i,t)
    y = zeros(T,1);
    p = zeros(N,T);
    x = zeros(N,T);
    w = zeros(N,T);
    v = zeros(N,T);
    W = zeros(N,T);
    U = zeros(N,T);
    nu = zeros(N,T);
    mu = zeros(N,T);
    for t=1:T
        y(t) = results.x(t);
        for i =1:N
            p(i,t) = results.x(T+(t-1)*N+i);
            x(i,t) = results.x(T+3*N*T+(t-1)*N+i);
            w(i,t) = results.x(T+N*T+(t-1)*N+i);
            v(i,t) = results.x(T+2*N*T+(t-1)*N+i);
            W(i,t) = results.x(T+4*N*T+(t-1)*N+i);
            U(i,t) = results.x(T+4*N*T+(t-1)*N+i);
            nu(i,t) = results.x(T+6*N*T+(t-1)*N+i);
            mu(i,t) = results.x(T+7*N*T+(t-1)*N+i);
        end
    end  
    obj = results.objval;
    
end