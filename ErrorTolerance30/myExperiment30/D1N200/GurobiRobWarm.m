function [obj,objBound,myGap]=GurobiRobWarm(myData,N,T,Zeta,solvec)
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
    NumOfInputs = 0; 
    %Objective constrints, Capacity constraints and ramping constraints
    NumOfRow=NumOfRow+2*((N+1)*T+N*T+T+T);
    NumOfRow = NumOfRow+2*N*T+(T-1)*N +N*T;
    NumOfInputs = NumOfInputs+2*(2*T*(N+1)+T*(N+1)*N+2*(N+1)*T+T*N*(N-1)/2+T*N+T*(N+1))+2*N*T+2*N*T+...
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
    InCount = 0;
    %Y(t),p(i,t), w(i,t), v(i,t),  x(i,t),W(i,t), U(i,t), nu(i,t),mu(i,t)
    %Objective Constraint 1(1)
    for t = 1:T
        for k = 1:N+1
            ii=ii+1;
            InCount = InCount+1;
            myRow(InCount) = ii;
            myCol(InCount) = t; 
            myInt(InCount) = 1; 
            %ConstraintBody(ii,t) = 1;
            
            InCount = InCount+1;
            myRow(InCount) = ii;
            myCol(InCount) = T+4*N*T+(t-1)*(N+1)+k; 
            myInt(InCount) = -1; 
            %ConstraintBody(ii,T+4*N*T+(t-1)*(N+1)+k) = -1;
            
            
            myRow(InCount+1:InCount+N) = ii;
            myCol(InCount+1:InCount+N) = T+(t-1)*N+(1:N); 
            myInt(InCount+1:InCount+N) = -myData.('HoldingCost')(t)+myData.('CommitmentCost')(1:N); 
            InCount = InCount+N;
            %ConstraintBody(ii,T+(t-1)*N+(1:N)) = -myData.('HoldingCost')(t)+myData.('CommitmentCost')(1:N);
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
            
            InCount = InCount+1;
            myRow(InCount) = ii;
            myCol(InCount) = T+4*N*T+(t-1)*(N+1)+k; 
            myInt(InCount) = 1; 
            %ConstraintBody(ii,T+4*N*T+(t-1)*(N+1)+k) = 1;
            
            InCount = InCount+1;
            myRow(InCount) = ii;
            myCol(InCount) = T+4*N*T+2*(N+1)*T+(t-1)*(N+1)+k; 
            myInt(InCount) = -Mt;  
            %ConstraintBody(ii,T+4*N*T+2*(N+1)*T+(t-1)*(N+1)+k) = -Mt;
            
            
            myRow(InCount+1:InCount+k-1) = ii;
            myCol(InCount+1:InCount+k-1) = T+(t-1)*N+(phi(t,1:k-1)); 
            myInt(InCount+1:InCount+k-1) = (-SortedH(t,1:k-1) + SortedH(t,k))...
                    .*gamH(t,phi(t,1:k-1));  
            InCount = InCount+k-1;
            %ConstraintBody(ii,T+(t-1)*N+(phi(t,1:k-1))) = (-SortedH(t,1:k-1) + SortedH(t,k))...
            %        .*gamH(t,phi(t,1:k-1));
%             for i = 1:k-1
%                 ConstraintBody(ii,T+(t-1)*N+phi(t,i)) = (-SortedH(t,i) + SortedH(t,k))...
%                     *gamH(t,phi(t,i));
%             end
        end
        ii = ii+1;
        InCount = InCount+1;
        myRow(InCount) = ii;
        myCol(InCount) = T+4*N*T+t*(N+1); 
        myInt(InCount) = 1;  
        %ConstraintBody(ii,T+4*N*T+t*(N+1)) = 1; 
        
        InCount = InCount+1;
        myRow(InCount) = ii;
        myCol(InCount) = T+4*N*T+2*(N+1)*T+t*(N+1); 
        myInt(InCount) = -Mt;  
        %ConstraintBody(ii,T+4*N*T+2*(N+1)*T+t*(N+1)) = -Mt;
        
        myRow(InCount+1:InCount+N) = ii;
        myCol(InCount+1:InCount+N) = T+(t-1)*N+phi(t,1:N); 
        myInt(InCount+1:InCount+N) = -SortedH(t,1:N).*gamH(t,phi(t,1:N));  
        InCount = InCount+N;
        %ConstraintBody(ii,T+(t-1)*N+phi(t,1:N)) =  -SortedH(t,1:N).*gamH(t,phi(t,1:N));
%         for i = 1:N
%             ConstraintBody(ii,T+(t-1)*N+phi(t,i)) =  -SortedH(t,i)*gamH(t,phi(t,i));
%         end
    end
    %Objective Constraint 1(3)
    for t = 1:T
        ii=ii+1;
        myRow(InCount+1:InCount+N+1) = ii;
        myCol(InCount+1:InCount+N+1) = T+4*N*T+2*(N+1)*T+(t-1)*(N+1)+(1:N+1); 
        myInt(InCount+1:InCount+N+1) = 1;  
        InCount = InCount+N+1;
        %ConstraintBody(ii,T+4*N*T+2*(N+1)*T+(t-1)*(N+1)+(1:N+1)) = 1;
%         for i = 1:N+1
%             ConstraintBody(ii,T+4*N*T+2*(N+1)*T+(t-1)*(N+1)+i) = 1; 
%         end
    end
    
    %Y(t),p(i,t), w(i,t), v(i,t),  x(i,t),W(i,t), U(i,t), nu(i,t),mu(i,t)
    %Objective Constraint 2(1)
    for t = 1:T
        for k = 1:N+1
            ii=ii+1;
            InCount = InCount+1;
            myRow(InCount) = ii;
            myCol(InCount) = t; 
            myInt(InCount) = 1;  
            %ConstraintBody(ii,t) = 1;
            
            InCount = InCount+1;
            myRow(InCount) = ii;
            myCol(InCount) = T+4*N*T+(N+1)*T+(t-1)*(N+1)+k; 
            myInt(InCount) = -1;  
            %ConstraintBody(ii,T+4*N*T+(N+1)*T+(t-1)*(N+1)+k) = -1;
            
            myRow(InCount+1:InCount+N) = ii;
            myCol(InCount+1:InCount+N) = T+(t-1)*N+(1:N); 
            myInt(InCount+1:InCount+N) = myData.('ShortageCost')(t)+myData.('CommitmentCost')(1:N);  
            InCount = InCount+N;
            %ConstraintBody(ii,T+(t-1)*N+(1:N)) = myData.('ShortageCost')(t)+myData.('CommitmentCost')(1:N);
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
            InCount = InCount+1;
            myRow(InCount) = ii;
            myCol(InCount) = T+4*N*T+(N+1)*T+(t-1)*(N+1)+k; 
            myInt(InCount) = 1;  
            %ConstraintBody(ii,T+4*N*T+(N+1)*T+(t-1)*(N+1)+k) = 1;
            
            InCount = InCount+1;
            myRow(InCount) = ii;
            myCol(InCount) = T+4*N*T+3*(N+1)*T+(t-1)*(N+1)+k; 
            myInt(InCount) = -Mt; 
            %ConstraintBody(ii,T+4*N*T+3*(N+1)*T+(t-1)*(N+1)+k) = -Mt;
            
            
            myRow(InCount+1:InCount+k-1) = ii;
            myCol(InCount+1:InCount+k-1) = T+(t-1)*N+psi(t,1:k-1); 
            myInt(InCount+1:InCount+k-1) = (-SortedS(t,1:k-1) + SortedS(t,k))...
                    .*etaS(t,psi(t,1:k-1)); 
            InCount = InCount+k-1;
            %ConstraintBody(ii,T+(t-1)*N+psi(t,1:k-1)) = (-SortedS(t,1:k-1) + SortedS(t,k))...
                    %.*etaS(t,psi(t,1:k-1));
%             for i = 1:k-1
%                 ConstraintBody(ii,T+(t-1)*N+psi(t,i)) = (-SortedS(t,i) + SortedS(t,k))...
%                     *etaS(t,psi(t,i));
%             end
        end
        ii = ii+1;
        InCount = InCount+1;
        myRow(InCount) = ii;
        myCol(InCount) = T+4*N*T+(N+1)*T+t*(N+1); 
        myInt(InCount) = 1;  
        %ConstraintBody(ii,T+4*N*T+(N+1)*T+t*(N+1)) = 1; 
        
        InCount = InCount+1;
        myRow(InCount) = ii;
        myCol(InCount) = T+4*N*T+3*(N+1)*T+t*(N+1); 
        myInt(InCount) = -Mt;  
        %ConstraintBody(ii,T+4*N*T+3*(N+1)*T+t*(N+1)) = -Mt;
        
        myRow(InCount+1:InCount+N) = ii;
        myCol(InCount+1:InCount+N) = T+(t-1)*N+psi(t,1:N); 
        myInt(InCount+1:InCount+N) =  -SortedS(t,1:N).*etaS(t,psi(t,1:N));
        InCount = InCount+N;
        %ConstraintBody(ii,T+(t-1)*N+psi(t,1:N)) = -SortedS(t,1:N).*etaS(t,psi(t,1:N));
%         for i = 1:N
%             ConstraintBody(ii,T+(t-1)*N+psi(t,i)) = -SortedS(t,i)*etaS(t,psi(t,i));
%         end
    end

    %Objective Constraint 2(3)
    %Y(t),p(i,t), w(i,t), v(i,t),  x(i,t),W(i,t), U(i,t), nu(i,t),mu(i,t)
    for t = 1:T
        ii=ii+1;
        myRow(InCount+1:InCount+N+1) = ii;
        myCol(InCount+1:InCount+N+1) = T+4*N*T+3*(N+1)*T+(t-1)*(N+1)+(1:N+1); 
        myInt(InCount+1:InCount+N+1) = 1;  
        InCount = InCount+N+1;
        %ConstraintBody(ii,T+4*N*T+3*(N+1)*T+(t-1)*(N+1)+(1:N+1)) = 1;
%         for i = 1:N+1
%             ConstraintBody(ii,T+4*N*T+3*(N+1)*T+(t-1)*(N+1)+i) = 1;
%         end
    end
     %Y(t), p(i,t), w(i,t), v(i,t),  x(i,t),W(i,t), U(i,t), nu(i,t),mu(i,t)
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
        i
        ConstraintBody(myRow(i),myCol(i)) = myInt(i);
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
    params.TimeLimit = 1000;
    results = gurobi(model,params);
    myGap = abs((results.objval - results.objbound)/results.objval);

    
    objBound = results.objbound;
    obj = results.objval;
    
end