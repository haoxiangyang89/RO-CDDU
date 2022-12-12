function [obj,p,detTime,solvec]=GurobiDet(myData,N,T)
    model.modelName = 'myDRModel';
    model.modelsense = 'Min';
    
    M=max(myData.('pmax'));
    
    %Sequence of decision varialbe in matrix 
    %Y(t),p(i,t), w(i,t), v(i,t),  x(i,t)
    %This step constructs the model objective
    myObjArr = [ones(T,1);zeros(4*N*T,1);]; 
    model.obj  = myObjArr;
    
    %Y(t),p(i,t), w(i,t), v(i,t),  x(i,t)
    %Set the lower bound for the decision variable
    DVLowerBound = [-inf*ones(T,1);zeros(4*N*T,1)];
    model.lb    = DVLowerBound;
    %Set the upper bound for the decision variables
    DVUpperBound = [inf*ones(T+N*T,1);ones(3*N*T,1)];
    model.ub    = DVUpperBound;

    %Set the datatype for the decision variables
    Datatype = [repmat('C',[1,T+N*T]),repmat('B',[1,3*N*T])];
    model.vtype = Datatype;

    %Count the total number of constraints
    NumOfCol = T+4*N*T;
    NumOfRow=0;
    NumOfInputs = 0; 
    %Objective constrints, Capacity constraints and ramping constraints
    NumOfRow=NumOfRow+T*2+2*N*T+(T-1)*N +N*T;
    NumOfInputs = NumOfInputs+T*(N+1)+T*(N+1)+2*N*T+2*N*T+...
        3*(T-1)*(N)+3*(T-1)*(N)+2*N; 
    
   NumOfInputs= NumOfInputs +3*(T-1)*(N)+2*N;
    %Upupmin constraints
    NumOfRow = NumOfRow+N*T;
    for i=1:N
        for t=1:T
            for tau=t:min(t+myData.('Upupmin')(i)-1,T)
                NumOfRow=NumOfRow+1;
                NumOfInputs = NumOfInputs+3;
            end
        end
    end
    NumOfInputs = NumOfInputs-N;
    
    NumOfInputs= NumOfInputs +3*(T-1)*(N)+2*N;
    %Downdownmin constraints
    NumOfRow = NumOfRow+N*T;
    for i=1:N
        for t=1:T
            for tau=t:min(t+myData.('Downdownmin')(i)-1,T)
                NumOfRow=NumOfRow+1;
                NumOfInputs = NumOfInputs+3;
            end
        end
    end
    NumOfInputs = NumOfInputs-N;
    
    %Set the main constraint matrix.
    ConstraintBody=sparse(NumOfRow,NumOfCol);   
    
    myRow = zeros(NumOfInputs,1);
    myCol = zeros(NumOfInputs,1);
    myInt = zeros(NumOfInputs,1);
    %Row index
    ii = 0;
    %Objective Constraint 1
    InCount = 0;
    for t = 1:T
        ii=ii+1;
        InCount = InCount+1;
        myRow(InCount) = ii;
        myCol(InCount) = t; 
        myInt(InCount) = 1; 
        %myA(ii,t) = 1;
        for i = 1:N 
            InCount = InCount+1;
            myRow(InCount) = ii;
            myCol(InCount) = T+(t-1)*N+i; 
            myInt(InCount) = -myData.('HoldingCost')(t)+myData.('CommitmentCost')(i); 
            %myA(ii,T+(t-1)*N+i) = -myData.('HoldingCost')(t)+myData.('CommitmentCost')(i);
        end
    end
    %Objective Constraint 2 
    for t = 1:T
        ii=ii+1;
        InCount = InCount+1;
        myRow(InCount) = ii;
        myCol(InCount) = t; 
        myInt(InCount) = 1; 
        %myA(ii,t) = 1;
        for i = 1:N 
          InCount = InCount+1;
          myRow(InCount) = ii;
          myCol(InCount) = T+(t-1)*N+i; 
          myInt(InCount) = myData.('ShortageCost')(t)+myData.('CommitmentCost')(i);           
           %myA(ii,T+(t-1)*N+i) = myData.('ShortageCost')(t)+myData.('CommitmentCost')(i);
        end
    end
     %Y(t), p(i,t), w(i,t), v(i,t),  x(i,t)
    %pmin
    for i =1:N
        for t=1:T
            ii=ii+1;
            InCount = InCount+1;
            myRow(InCount) = ii;
            myCol(InCount) = T+(t-1)*N+i; 
            myInt(InCount) = 1; 
            %myA(ii,T+(t-1)*N+i)=1;
            
            InCount = InCount+1;
            myRow(InCount) = ii;
            myCol(InCount) = T+3*N*T+(t-1)*N+i; 
            myInt(InCount) = -myData.('pmin')(i); 
            %myA(ii,T+3*N*T+(t-1)*N+i)=-myData.('pmin')(i);
        end
    end

    %Y(t), p(i,t), w(i,t), v(i,t),  x(i,t)
    %pmax
    for i =1:N
        for t=1:T
            ii=ii+1;
            InCount = InCount+1;
            myRow(InCount) = ii;
            myCol(InCount) = T+(t-1)*N+i; 
            myInt(InCount) = -1; 
            %myA(ii,T+(t-1)*N+i)=-1;
            
            InCount = InCount+1;
            myRow(InCount) = ii;
            myCol(InCount) = T+3*N*T+(t-1)*N+i; 
            myInt(InCount) = myData.('pmax')(i); 
            %myA(ii,T+3*N*T+(t-1)*N+i)=myData.('pmax')(i);
        end
    end
    %Y(t), p(i,t), w(i,t), v(i,t),  x(i,t)
    %RD
    for i = 1:N
        for t= 1:T-1
            ii=ii+1;
            InCount = InCount+1;
            myRow(InCount) = ii;
            myCol(InCount) = T+t*N+i; 
            myInt(InCount) = 1; 
            %myA(ii,T+t*N+i)=1;
            
            InCount = InCount+1;
            myRow(InCount) = ii;
            myCol(InCount) = T+(t-1)*N+i; 
            myInt(InCount) = -1; 
            %myA(ii,T+(t-1)*N+i)=-1;
            
            InCount = InCount+1;
            myRow(InCount) = ii;
            myCol(InCount) =T+3*N*T+(t-1)*N+i; 
            myInt(InCount) = myData.('RD')(i); 
            %myA(ii,T+3*N*T+(t-1)*N+i) = myData.('RD')(i);
        end
    end
    %Y(t), p(i,t), w(i,t), v(i,t),  x(i,t)
    %RU
    for i = 1:N
        ii=ii+1;
        InCount = InCount+1;
        myRow(InCount) = ii;
        myCol(InCount) =T+ii; 
        myInt(InCount) = -1; 
        %myA(ii,T+i)=-1;
        
        InCount = InCount+1;
        myRow(InCount) = ii;
        myCol(InCount) =T+3*N*T+i; 
        myInt(InCount) = myData.('RU')(i); 
        %myA(ii,T+3*N*T+i)=myData.('RU')(i);
        for t=2:T
            ii=ii+1;
            InCount = InCount+1;
            myRow(InCount) = ii;
            myCol(InCount) =T+(t-1)*N+i; 
            myInt(InCount) = -1; 
            %myA(ii,T+(t-1)*N+i)=-1;
            
            InCount = InCount+1;
            myRow(InCount) = ii;
            myCol(InCount) =T+(t-2)*N+i; 
            myInt(InCount) = 1; 
            %myA(ii,T+(t-2)*N+i)=1;
            
            InCount = InCount+1;
            myRow(InCount) = ii;
            myCol(InCount) =T+3*N*T+(t-1)*N+i; 
            myInt(InCount) = myData.('RU')(i); 
            %myA(ii,T+3*N*T+(t-1)*N+i)=myData.('RU')(i);
        end
    end
    
    %Y(t), p(i,t), w(i,t), v(i,t),  x(i,t)
    %Upupmin
    for i= 1:N
        ii=ii+1;
        InCount = InCount+1;
        myRow(InCount) = ii;
        myCol(InCount) =T+i; 
        myInt(InCount) = -1; 
        %myA(ii,T+i) = -1;
        
        InCount = InCount+1;
        myRow(InCount) = ii;
        myCol(InCount) = T+N*T+i; 
        myInt(InCount) = M; 
        %myA(ii,T+N*T+i) = M;
        for t=2:T
            ii=ii+1;
            InCount = InCount+1;
            myRow(InCount) = ii;
            myCol(InCount) = T+(t-1)*N+i; 
            myInt(InCount) = -1; 
            %myA(ii,T+(t-1)*N+i)=-1;
            
            InCount = InCount+1;
            myRow(InCount) = ii;
            myCol(InCount) = T+(t-2)*N+i; 
            myInt(InCount) = 1; 
            %myA(ii,T+(t-2)*N+i)=1;
            
            InCount = InCount+1;
            myRow(InCount) = ii;
            myCol(InCount) = T+N*T+(t-1)*N+i; 
            myInt(InCount) = M; 
            %myA(ii,T+N*T+(t-1)*N+i)=M;
        end
    end
    for i=1:N
        ii=ii+1;
        InCount = InCount+1;
        myRow(InCount) = ii;
        myCol(InCount) = T+i; 
        myInt(InCount) = 1; 
        %myA(ii,T+i)=1;
        
        InCount = InCount+1;
        myRow(InCount) = ii;
        myCol(InCount) = T+N*T+i; 
        myInt(InCount) = -M; 
        %myA(ii,T+N*T+i)=-M;
        
        for tau = 2:min(myData.('Upupmin')(i),T)
            ii=ii+1;
            InCount = InCount+1;
            myRow(InCount) = ii;
            myCol(InCount) = T+(tau-1)*N+i; 
            myInt(InCount) = 1; 
            %myA(ii,T+(tau-1)*N+i)=1;
            
            InCount = InCount+1;
            myRow(InCount) = ii;
            myCol(InCount) = T+(tau-2)*N+i; 
            myInt(InCount) = -1; 
            %myA(ii,T+(tau-2)*N+i)=-1;
            
            InCount = InCount+1;
            myRow(InCount) = ii;
            myCol(InCount) = T+N*T+i; 
            myInt(InCount) = -M; 
            %myA(ii,T+N*T+i)=-M;
        end
        for t= 2: T 
            for tau = t:min(t+myData.('Upupmin')(i)-1,T)
                ii=ii+1;
                InCount = InCount+1;
                myRow(InCount) = ii;
                myCol(InCount) = T+(tau-1)*N+i; 
                myInt(InCount) = 1; 
                %myA(ii,T+(tau-1)*N+i)=1;
                
                InCount = InCount+1;
                myRow(InCount) = ii;
                myCol(InCount) =T+(tau-2)*N+i; 
                myInt(InCount) = -1; 
                %myA(ii,T+(tau-2)*N+i)=-1;
                
                InCount = InCount+1;
                myRow(InCount) = ii;
                myCol(InCount) =T+N*T+(t-1)*N+i; 
                myInt(InCount) = -M; 
                %myA(ii,T+N*T+(t-1)*N+i)=-M;
            end
        end
    end

    %Y(t), p(i,t), w(i,t), v(i,t),  x(i,t)
    %Downdownmin
    for i= 1:N
        ii=ii+1;
        InCount = InCount+1;
        myRow(InCount) = ii;
        myCol(InCount) =T+i;
        myInt(InCount) = 1; 
       % myA(ii,T+i) = 1;
        
        InCount = InCount+1;
        myRow(InCount) = ii;
        myCol(InCount) =T+2*N*T+i;
        myInt(InCount) = M; 
        %myA(ii,T+2*N*T+i) = M;
        for t=2:T
            ii=ii+1;
            InCount = InCount+1;
            myRow(InCount) = ii;
            myCol(InCount) =T+(t-1)*N+i;
            myInt(InCount) = 1; 
            %myA(ii,T+(t-1)*N+i)=1;
            
            InCount = InCount+1;
            myRow(InCount) = ii;
            myCol(InCount) =T+(t-2)*N+i;
            myInt(InCount) = -1; 
            %myA(ii,T+(t-2)*N+i)=-1;
            
            InCount = InCount+1;
            myRow(InCount) = ii;
            myCol(InCount) =T+2*N*T+(t-1)*N+i;
            myInt(InCount) = M; 
            %myA(ii,T+2*N*T+(t-1)*N+i)=M;
        end
    end
    for i =1:N
        ii=ii+1;
        InCount = InCount+1;
        myRow(InCount) = ii;
        myCol(InCount) =T+i;
        myInt(InCount) = -1; 
        %myA(ii,T+i)=-1;
        
        InCount = InCount+1;
        myRow(InCount) = ii;
        myCol(InCount) =T+2*N*T+i;
        myInt(InCount) = -M; 
        %myA(ii,T+2*N*T+i)=-M;
        for tau = 2:min(myData.('Downdownmin')(i),T)
            ii=ii+1;
            InCount = InCount+1;
            myRow(InCount) = ii;
            myCol(InCount) =T+(tau-1)*N+i;
            myInt(InCount) = -1; 
           % myA(ii,T+(tau-1)*N+i)=-1;
            
            InCount = InCount+1;
            myRow(InCount) = ii;
            myCol(InCount) =T+(tau-2)*N+i;
            myInt(InCount) = 1; 
            %myA(ii,T+(tau-2)*N+i)=1;
            
            InCount = InCount+1;
            myRow(InCount) = ii;
            myCol(InCount) =T+(tau-2)*N+i;
            myInt(InCount) = 1;
            %myA(ii,T+(tau-2)*N+i)= 1;
        end
        for t = 2:T
            for tau = t:min(t+myData.('Downdownmin')(i)-1,T)
                ii=ii+1;
                InCount = InCount+1;
                myRow(InCount) = ii;
                myCol(InCount) =T+(tau-1)*N+i;
                myInt(InCount) = -1;
                %myA(ii,T+(tau-1)*N+i)=-1;
                
                InCount = InCount+1;
                myRow(InCount) = ii;
                myCol(InCount) =T+(tau-2)*N+i;
                myInt(InCount) = 1;
                %myA(ii,T+(tau-2)*N+i)=1;
                
                InCount = InCount+1;
                myRow(InCount) = ii;
                myCol(InCount) =T+2*N*T+(t-1)*N+i;
                myInt(InCount) = -M;
                %myA(ii,T+2*N*T+(t-1)*N+i)=-M;
            end
        end
    end
    %ConstraintBody(1:end,1:end) = myA;
    for i = 1:NumOfInputs
        ConstraintBody(myRow(i),myCol(i)) = myInt(i);
    end
    model.A =  ConstraintBody;
    
    
    %Constraintbody's bounds on A
    Alower=zeros(NumOfRow,1);
    
    ii=0;
    %Objective Constraint 1 
    for t = 1:T
        ii=ii+1;
        Alower(ii,1)=-myData.('HoldingCost')(t)*myData.('D')(t); 
    end
   
    %Objective Constraint 2
    for t =1:T
        ii=ii+1;
        Alower(ii,1)=myData.('ShortageCost')(t)*myData.('D')(t);
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
    params.MIPGap = 0.01;
    params.TimeLimit = 1000;
 
    tic
    results = gurobi(model,params);
    detTime=toc;
    
    %Y(t), p(i,t), w(i,t), v(i,t),  x(i,t)
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
    solvec = results.x;
end