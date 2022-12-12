function [obj,p,solvec,ImprUppTime]=GurobiMcCormick(myData,N,T,Zeta)
%try
    %Sorted value  
    %Sorted value  
    %Construct the parameter for later model
    SortedH = zeros(T,N+1);
    SortedS = zeros(T,N+1);
    %Sorted index
    phi = zeros(T,N);
    psi = zeros(T,N);
    gamH = zeros(T,N);
    etaS = zeros(T,N);
    
    %Upperbound for pi and lambda
    piupper = zeros(T*N,1);
    lambdaupper = zeros(T*N,1);
    for t = 1:T
        myH = abs(myData.('HoldingCost')(t) - myData.('CommitmentCost'));
        [SortedH(t,1:N),phi(t,:)] = sort(myH,'descend');
        piupper((t-1)*N+1:t*N)= myH';
        myS = abs(myData.('ShortageCost')(t)+myData.('CommitmentCost'));
        [SortedS(t,1:N),psi(t,:)] = sort(myS,'descend');
        lambdaupper((t-1)*N+1:t*N) = myS';
        
        ck = (myData.('HoldingCost')(t)<= myData.('CommitmentCost'));
        gamH(t,:) = abs(myData.('alpha')).*ck + abs(myData.('beta')).*(1-ck);
        ck = (-myData.('ShortageCost')(t)<=myData.('CommitmentCost'));
        etaS(t,:) = abs(myData.('alpha')).*ck + abs(myData.('beta')).*(1-ck);
        
    end
    
    
    model.modelName = 'myDRModelMcCormick';
    model.modelsense = 'Min';

    M =max(myData.('pmax'));
    
    %Sequence of decision varialbe in matrix 
    %Y(t), p(i,t), w(i,t), v(i,t),  x(i,t), pi(i,t), lambda(i,t),
    %mu(t),theta(t), pi*p(i,t),eta*p(i,t) 
    %This step constructs the model objective
    myObjArr = [ones(T,1);zeros(4*N*T+4*N*T+2*T,1)]; 
    model.obj   = myObjArr;

    %Y(t), p(i,t), w(i,t), v(i,t),  x(i,t), pi(i,t), lambda(i,t),
    %mu(t),theta(t), pi*p(i,t),eta*p(i,t) 
    %Set the lower bound for the decision variable
    DVLowerBound = [-inf*ones(T,1);zeros(4*N*T,1);zeros(4*N*T+2*T,1)];
    model.lb    = DVLowerBound;
    %Set the upper bound for the decision variables
    DVUpperBound = [inf*ones(T+N*T,1);ones(3*N*T,1);piupper;lambdaupper;SortedH(:,1);SortedS(:,1);inf*ones(2*N*T,1)];
    model.ub    = DVUpperBound;

    %Set the datatype for the decision variables
    Datatype = [repmat('C',[1,T+N*T]),repmat('B',[1,3*N*T]),repmat('C',[1,4*N*T+2*T])];
    model.vtype = Datatype;

    %Set the main constraint matrix.
    NumOfRow = 0;
    NumOfCol = T+4*N*T+4*N*T+2*T;
    NumOfInputs = 0; 


    %Objective contraints and Capacity constraints and ramping constraints 
    NumOfRow=NumOfRow+T*2+2*N*T+(T-1)*N+N*T;
    NumOfInputs = NumOfInputs+2*T*(2*N+2)+2*N*T+2*N*T+...
        3*(T-1)*(N)+3*(T-1)*(N)+2*N;

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
    
    %Dual problem constraints 
    NumOfRow = NumOfRow+2*N*T;
    NumOfInputs = NumOfInputs+2*N*T*2;
    
    %McCormick Constraints
    NumOfRow = NumOfRow + 6*N*T;
    NumOfInputs = NumOfInputs + 2*(2*2*N*T+3*N*T);
    
    %Set the main constraint matrix.
    ConstraintBody=sparse(NumOfRow,NumOfCol);   
    myRow = zeros(NumOfInputs,1);
    myCol = zeros(NumOfInputs,1);
    myInt = zeros(NumOfInputs,1);

    ii=0;
    InCount = 0;
    %Y(t), p(i,t), w(i,t), v(i,t),  x(i,t), pi(i,t), lambda(i,t),
    %mu(t),theta(t), pi*p(i,t),eta*p(i,t) 
    %Objective Constraint 1 
    for t = 1:T
        ii=ii+1;
        InCount = InCount+1;
        myRow(InCount) = ii;
        myCol(InCount) = t; 
        myInt(InCount) = 1; 
        %ConstraintBody(ii,t) = 1;
        
        myRow(InCount+1:InCount+N) = ii;
        myCol(InCount+1:InCount+N) =T+(t-1)*N+(1:N); 
        myInt(InCount+1:InCount+N) = -myData.('HoldingCost')(t)+myData.('CommitmentCost')(1:N); 
        InCount = InCount+N;
        
        InCount = InCount+1;
        myRow(InCount) =ii;
        myCol(InCount) = T+6*N*T+t;
        myInt(InCount) = -Zeta;
        
        myRow(InCount+1:InCount+N)= ii; 
        myCol(InCount+1:InCount+N) = 3*T+6*N*T+(t-1)*N+(1:N); 
        myInt(InCount+1:InCount+N) = -(gamH(t,1:N))';
        InCount = InCount +N;
%         for i = 1:N 
%             ConstraintBody(ii,T+(t-1)*N+i) = -gamH(i,t)*Pi(i,t)-myData.('HoldingCost')(t)+myData.('CommitmentCost')(i);
%         end
    end

    %Y(t), p(i,t), w(i,t), v(i,t),  x(i,t), pi(i,t), lambda(i,t),
    %mu(t),theta(t), pi*p(i,t),eta*p(i,t) 
    %Objective Constraint 2
    for t= 1:T
        ii=ii+1;
        InCount = InCount+1;
        myRow(InCount) = ii;
        myCol(InCount) = t; 
        myInt(InCount) = 1; 
        %ConstraintBody(ii,t) = 1;
        
        myRow(InCount+1:InCount+N) = ii;
        myCol(InCount+1:InCount+N) = T+(t-1)*N+(1:N); 
        myInt(InCount+1:InCount+N) = myData.('ShortageCost')(t)+myData.('CommitmentCost')(1:N); 
        InCount = InCount+N;
        
        
        InCount = InCount+1;
        myRow(InCount) =ii;
        myCol(InCount) = T+6*N*T+T+t;
        myInt(InCount) = -Zeta;
        
        myRow(InCount+1:InCount+N)= ii; 
        myCol(InCount+1:InCount+N) =3*T+7*N*T+(t-1)*N+(1:N); 
        myInt(InCount+1:InCount+N) = -(etaS(t,1:N))';
        InCount = InCount +N;
%         for i=1:N
%             ConstraintBody(ii,T+(t-1)*N+i) = -etaS(i,t)*Lam(i,t)+myData.('ShortageCost')(t)+myData.('CommitmentCost')(i);
%         end
    end

    %Y(t), p(i,t), w(i,t), v(i,t),  x(i,t), pi(i,t), lambda(i,t),
    %mu(t),theta(t), pi*p(i,t),eta*p(i,t) 
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

    %Y(t), p(i,t), w(i,t), v(i,t),  x(i,t), pi(i,t), lambda(i,t),
    %mu(t),theta(t), pi*p(i,t),eta*p(i,t) 
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
    %Y(t), p(i,t), w(i,t), v(i,t),  x(i,t), pi(i,t), lambda(i,t),
    %mu(t),theta(t), pi*p(i,t),eta*p(i,t) 
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
    %Y(t), p(i,t), w(i,t), v(i,t),  x(i,t), pi(i,t), lambda(i,t),
    %mu(t),theta(t), pi*p(i,t),eta*p(i,t) 
    %RU
    for i = 1:N
        ii=ii+1;
        InCount = InCount+1;
        myRow(InCount) = ii;
        myCol(InCount) =T+i; 
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
    
    %Y(t), p(i,t), w(i,t), v(i,t),  x(i,t), pi(i,t), lambda(i,t),
    %mu(t),theta(t), pi*p(i,t),eta*p(i,t) 
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

    %Y(t), p(i,t), w(i,t), v(i,t),  x(i,t), pi(i,t), lambda(i,t),
    %mu(t),theta(t), pi*p(i,t),eta*p(i,t) 
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
    
    %Dual problem constraints 1
    %Y(t), p(i,t), w(i,t), v(i,t),  x(i,t), pi(i,t), lambda(i,t),
    %mu(t),theta(t), pi*p(i,t),eta*p(i,t) 
    for i = 1:N
        for t = 1:T
            ii=ii+1;
            InCount = InCount+1;
            myRow(InCount) = ii;
            myCol(InCount) = T+4*N*T+(t-1)*N+i;
            myInt(InCount) = 1;
            
            InCount = InCount+1;
            myRow(InCount) = ii;
            myCol(InCount) = T+6*N*T+t;
            myInt(InCount) = 1; 
        end
    end
    %Dual problem constraints 2
    %Y(t), p(i,t), w(i,t), v(i,t),  x(i,t), pi(i,t), lambda(i,t),
    %mu(t),theta(t), pi*p(i,t),eta*p(i,t) 
    for i = 1:N
        for t = 1:T
            ii=ii+1;
            InCount = InCount+1;
            myRow(InCount) = ii;
            myCol(InCount) = T+5*N*T+(t-1)*N+i;
            myInt(InCount) = 1;
            
            InCount = InCount+1;
            myRow(InCount) = ii;
            myCol(InCount) = T+6*N*T+T+t;
            myInt(InCount) = 1; 
        end
    end
    
    %McCormick Constraints 1(1)
    %Y(t), p(i,t), w(i,t), v(i,t),  x(i,t), pi(i,t), lambda(i,t),
    %mu(t),theta(t), pi*p(i,t),eta*p(i,t) 
    for i = 1:N
        for t=  1:T
            ii=ii+1;
            InCount = InCount+1;
            myRow(InCount) = ii;
            myCol(InCount) = 3*T+6*N*T+(t-1)*N+i;
            myInt(InCount) = -1; 
            
            InCount = InCount+1;
            myRow(InCount) = ii;
            myCol(InCount) = T+(t-1)*N+i;
            myInt(InCount) = piupper((t-1)*N+i); 
        end
    end
    %McCormick Constraints 1(2)
    %Y(t), p(i,t), w(i,t), v(i,t),  x(i,t), pi(i,t), lambda(i,t),
    %mu(t),theta(t), pi*p(i,t),eta*p(i,t) 
    for i = 1:N
        for t=  1:T
            ii=ii+1;
            InCount = InCount+1;
            myRow(InCount) = ii;
            myCol(InCount) = 3*T+6*N*T+(t-1)*N+i;
            myInt(InCount) = -1; 
            
            InCount = InCount+1;
            myRow(InCount) = ii;
            myCol(InCount) = T+4*N*T+(t-1)*N+i;
            myInt(InCount) = myData.('pmax')(i); 
        end
    end
    %McCormick Constraints 1(3)
    %Y(t), p(i,t), w(i,t), v(i,t),  x(i,t), pi(i,t), lambda(i,t),
    %mu(t),theta(t), pi*p(i,t),eta*p(i,t) 
    for i = 1:N
        for t=  1:T
            ii=ii+1;
            InCount = InCount+1;
            myRow(InCount) = ii;
            myCol(InCount) = 3*T+6*N*T+(t-1)*N+i;
            myInt(InCount) = 1; 
            
            InCount = InCount+1;
            myRow(InCount) = ii;
            myCol(InCount) = T+4*N*T+(t-1)*N+i;
            myInt(InCount) = -myData.('pmax')(i); 
            
            InCount = InCount+1;
            myRow(InCount) = ii;
            myCol(InCount) = T+(t-1)*N+i;
            myInt(InCount) = -piupper((t-1)*N+i); 
        end
    end
    %McCormick Constraints 2(1)
    %Y(t), p(i,t), w(i,t), v(i,t),  x(i,t), pi(i,t), lambda(i,t),
    %mu(t),theta(t), pi*p(i,t),eta*p(i,t) 
    for i = 1:N
        for t=  1:T
            ii=ii+1;
            InCount = InCount+1;
            myRow(InCount) = ii;
            myCol(InCount) = 3*T+7*N*T+(t-1)*N+i;
            myInt(InCount) = -1; 
            
            InCount = InCount+1;
            myRow(InCount) = ii;
            myCol(InCount) = T+(t-1)*N+i;
            myInt(InCount) = lambdaupper((t-1)*N+i); 
        end
    end
    %McCormick Constraints 2(2)
    %Y(t), p(i,t), w(i,t), v(i,t),  x(i,t), pi(i,t), lambda(i,t),
    %mu(t),theta(t), pi*p(i,t),eta*p(i,t) 
    for i = 1:N
        for t=  1:T
            ii=ii+1;
            InCount = InCount+1;
            myRow(InCount) = ii;
            myCol(InCount) = 3*T+7*N*T+(t-1)*N+i;
            myInt(InCount) = -1; 
            
            InCount = InCount+1;
            myRow(InCount) = ii;
            myCol(InCount) = T+5*N*T+(t-1)*N+i;
            myInt(InCount) = myData.('pmax')(i); 
        end
    end
    %McCormick Constraints 2(3)
    %Y(t), p(i,t), w(i,t), v(i,t),  x(i,t), pi(i,t), lambda(i,t),
    %mu(t),theta(t), pi*p(i,t),eta*p(i,t) 
    for i = 1:N
        for t=  1:T
            ii=ii+1;
            InCount = InCount+1;
            myRow(InCount) = ii;
            myCol(InCount) = 3*T+7*N*T+(t-1)*N+i;
            myInt(InCount) = 1; 
            
            InCount = InCount+1;
            myRow(InCount) = ii;
            myCol(InCount) = T+5*N*T+(t-1)*N+i;
            myInt(InCount) = -myData.('pmax')(i); 
            
            InCount = InCount+1;
            myRow(InCount) = ii;
            myCol(InCount) = T+(t-1)*N+i;
            myInt(InCount) = -lambdaupper((t-1)*N+i); 
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
    
   %Dual problem constraints 1
    %Y(t), p(i,t), w(i,t), v(i,t),  x(i,t), pi(i,t), lambda(i,t),
    %mu(t),theta(t), pi*p(i,t),eta*p(i,t) 
    for i = 1:N
        for t = 1:T
            ii=ii+1;
            Alower(ii,1)=piupper((t-1)*N+i);
        end
    end
    %Dual problem constraints 2
    %Y(t), p(i,t), w(i,t), v(i,t),  x(i,t), pi(i,t), lambda(i,t),
    %mu(t),theta(t), pi*p(i,t),eta*p(i,t) 
    for i = 1:N
        for t = 1:T
            ii=ii+1;
            Alower(ii,1)=lambdaupper((t-1)*N+i);
        end
    end
    
    %McCormick Constraints 1(1)
    %Y(t), p(i,t), w(i,t), v(i,t),  x(i,t), pi(i,t), lambda(i,t),
    %mu(t),theta(t), pi*p(i,t),eta*p(i,t) 
    for i = 1:N
        for t=  1:T
            ii=ii+1;
            Alower(ii,1) = 0;
        end
    end
    %McCormick Constraints 1(2)
    %Y(t), p(i,t), w(i,t), v(i,t),  x(i,t), pi(i,t), lambda(i,t),
    %mu(t),theta(t), pi*p(i,t),eta*p(i,t) 
    for i = 1:N
        for t=  1:T
            ii=ii+1;
            Alower(ii,1) = 0;
        end
    end
    %McCormick Constraints 1(3)
    %Y(t), p(i,t), w(i,t), v(i,t),  x(i,t), pi(i,t), lambda(i,t),
    %mu(t),theta(t), pi*p(i,t),eta*p(i,t) 
    for i = 1:N
        for t=  1:T
            ii=ii+1;
            Alower(ii,1) =-piupper((t-1)*N+i)*myData.('pmax')(i);
        end
    end
    %McCormick Constraints 2(1)
    %Y(t), p(i,t), w(i,t), v(i,t),  x(i,t), pi(i,t), lambda(i,t),
    %mu(t),theta(t), pi*p(i,t),eta*p(i,t) 
    for i = 1:N
        for t=  1:T
            ii=ii+1;
            Alower(ii,1) = 0; 
        end
    end
    %McCormick Constraints 2(2)
    %Y(t), p(i,t), w(i,t), v(i,t),  x(i,t), pi(i,t), lambda(i,t),
    %mu(t),theta(t), pi*p(i,t),eta*p(i,t) 
    for i = 1:N
        for t=  1:T
            ii=ii+1;
            Alower(ii,1) = 0;
        end
    end
    %McCormick Constraints 2(3)
    %Y(t), p(i,t), w(i,t), v(i,t),  x(i,t), pi(i,t), lambda(i,t),
    %mu(t),theta(t), pi*p(i,t),eta*p(i,t) 
    for i = 1:N
        for t=  1:T
            ii=ii+1;
            Alower(ii,1) =-lambdaupper((t-1)*N+i)*myData.('pmax')(i);
        end
    end
    
    
    
    model.rhs = Alower;
    model.sense = repmat('>',NumOfRow,1);
    
    clear params;
    params.MIPGap = 0.005;
    params.TimeLimit = 10000;
 
    tic
    results = gurobi(model,params);
    ImprUppTime = toc;
    
    %Y(t), p(i,t), w(i,t), v(i,t),  x(i,t), pi(i,t), eta(i,t),
    %mu(t),theta(t), pi*p(i,t),eta*p(i,t) 
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
%catch err
    %throw(err);
%end
end