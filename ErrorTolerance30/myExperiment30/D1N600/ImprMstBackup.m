function [obj,p,solvec,ImprUppTime]=ImprMst(myData,N,T,Pi,Lam,Mu,Tht,Zeta)
try
    %Sorted value  
    gamH = zeros(N,T);
    etaS = zeros(N,T);
    for t = 1:T
        ck = (myData.('HoldingCost')(t)<= myData.('CommitmentCost'));
        gamH(:,t) = transpose(abs(myData.('alpha')).*ck + abs(myData.('beta')).*(1-ck));
        ck = (-myData.('ShortageCost')(t)<=myData.('CommitmentCost'));
        etaS(:,t) = transpose(abs(myData.('alpha')).*ck + abs(myData.('beta')).*(1-ck));
    end
    
    
    model.modelName = 'myDRModel';
    model.modelsense = 'Min';

    M =max(myData.('pmax'));;
    
    %Sequence of decision varialbe in matrix 
    %Y(t), p(i,t), w(i,t), v(i,t),  x(i,t)
    %This step constructs the model objective
    myObjArr = [ones(T,1);zeros(4*N*T,1)]; 
    model.obj   = myObjArr;

    %Y(t), p(i,t), w(i,t), v(i,t),  x(i,t)
    %Set the lower bound for the decision variable
    DVLowerBound = [-inf*ones(T,1);zeros(4*N*T,1)];
    model.lb    = DVLowerBound;
    %Set the upper bound for the decision variables
    DVUpperBound = [inf*ones(T+N*T,1);ones(3*N*T,1)];
    model.ub    = DVUpperBound;

    %Set the datatype for the decision variables
    Datatype = [repmat('C',[1,T+N*T]),repmat('B',[1,3*N*T])];
    model.vtype = Datatype;

    %Set the main constraint matrix.
    NumOfRow = 0;
    NumOfCol = T+4*N*T;

    %Objective contraints and Capacity constraints and ramping constraints 
    NumOfRow=T*2+2*N*T+(T-1)*N+N*T;

    %Upupmin constraints
    NumOfRow =NumOfRow + N*T;
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


    ii=0;
    %Y(t), p(i,t), w(i,t), v(i,t),  x(i,t)
    %Objective Constraint 1 
    for t = 1:T
        ii=ii+1;
        ConstraintBody(ii,t) = 1;
        for i = 1:N 
            ConstraintBody(ii,T+(t-1)*N+i) = -gamH(i,t)*Pi(i,t)-myData.('HoldingCost')(t)+myData.('CommitmentCost')(i);
        end
    end

    %Y(t), p(i,t), w(i,t), v(i,t),  x(i,t)
    %Objective Constraint 2
    for t= 1:T
        ii=ii+1;
        ConstraintBody(ii,t) = 1;
        for i=1:N
            ConstraintBody(ii,T+(t-1)*N+i) = -etaS(i,t)*Lam(i,t)+myData.('ShortageCost')(t)+myData.('CommitmentCost')(i);
        end
    end

     %Y(t), p(i,t), w(i,t), v(i,t),  x(i,t)
    %pmin
    for i =1:N
        for t=1:T
            ii=ii+1;
            ConstraintBody(ii,T+(t-1)*N+i)=1;
            ConstraintBody(ii,T+3*N*T+(t-1)*N+i)=-myData.('pmin')(i);
        end
    end

     %Y(t), p(i,t), w(i,t), v(i,t),  x(i,t)
    %pmax
    for i =1:N
        for t=1:T
            ii=ii+1;
            ConstraintBody(ii,T+(t-1)*N+i)=-1;
            ConstraintBody(ii,T+3*N*T+(t-1)*N+i)=myData.('pmax')(i);
        end
    end
        %Y(t), p(i,t), w(i,t), v(i,t),  x(i,t)
    %RD
    for i = 1:N
        for t= 1:T-1
            ii=ii+1;
            ConstraintBody(ii,T+t*N+i)=1;
            ConstraintBody(ii,T+(t-1)*N+i)=-1;
            ConstraintBody(ii,T+3*N*T+(t-1)*N+i) = myData.('RD')(i);
        end
    end
    %Y(t), p(i,t), w(i,t), v(i,t),  x(i,t)
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
    
    %Y(t), p(i,t), w(i,t), v(i,t),  x(i,t)
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

    %Y(t), p(i,t), w(i,t), v(i,t),  x(i,t)
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
    %Objective Constraint 1 
    for t = 1:T
        ii=ii+1;
        Alower(ii,1)=Zeta*Mu(t)-myData.('HoldingCost')(t)*myData.('D')(t); 
    end
   
    %Objective Constraint 2
    for t =1:T
        ii=ii+1;
        Alower(ii,1)=Zeta*Tht(t)+myData.('ShortageCost')(t)*myData.('D')(t);
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

    model.rhs = Alower;
    model.sense = repmat('>',NumOfRow,1);
    
    clear params;
    params.MIPGap = 0.005;
 
    tic
    results = gurobi(model,params);
    ImprUppTime = toc;
    
    %Y(t), p(i,t), w(i,t), v(i,t),  x(i,t)
    solvec = results.x;
    y = zeros(T,1);
    p = zeros(N,T);
    x = zeros(N,T);
    w = zeros(N,T);
    v = zeros(N,T);
    for t=1:T
        y(t) = results.x(t);
        for i =1:N
            p(i,t) = results.x(T+(t-1)*N+i);
            x(i,t) = results.x(T+3*N*T+(t-1)*N+i);
            w(i,t) = results.x(T+N*T+(t-1)*N+i);
            v(i,t) = results.x(T+2*N*T+(t-1)*N+i);     
        end
    end  
    obj = results.objval;
catch err
    throw(err);
end
end