%%% Designed, Zhang Qi,Dong Yingjie,Ye shan,Li Xu,He Dongcheng,Xiang Guoqi %%%
function[Best_score,Best_pos,TNTWCOA_curve,Trajectories,fitness_history,position_history]=TNTWCOA(SearchAgents,Max_iterations,lowerbound,upperbound,dimension,fitness)


Trajectories=zeros(SearchAgents,Max_iterations);
position_history=zeros(SearchAgents,Max_iterations,dimension);
fitness_history=zeros(SearchAgents,Max_iterations);


lowerbound=ones(1,dimension).*(lowerbound);                              % Lower limit for variables
upperbound=ones(1,dimension).*(upperbound);                              % Upper limit for variables

%% Chaotic mapping strategy for Algorithm initialization process %%
%The chaotic initialization function chaos(N,SearchAgents,dimension), where N is 1-10, represents tent mapping, Logistic mapping, Cubic mapping, chebyshev mapping, Piecewise mapping, sinusoidal mapping, Sine mapping, respectively. ICMIC mapping, Circle mapping,Bernoulli mapping.
X = repmat(lowerbound,SearchAgents,1)+chaos(1,SearchAgents,dimension).* repmat((upperbound-lowerbound),SearchAgents,1);%


for i =1:SearchAgents
    L=X(i,:);
    fit(i)=fitness(L);
end


for t=1:Max_iterations
     for ii=1:SearchAgents
        position_history(ii,t,:)=X(ii,:);
        Trajectories(:,t)=X(:,1);
        fitness_history(ii,t)=fit(1,ii);
    end


    %% update the best condidate solution%%
    [best , location]=min(fit);
    if t==1
        Xbest=X(location,:);                                           % Optimal location
        fbest=best;                                           % The optimization objective function
    elseif best<fbest
        fbest=best;
        Xbest=X(location,:);
    end

    %%Nonlinear Inertia Weight Factor for Hunting and attacking strategy on iguana%%
    
    w = (exp(t/Max_iterations)-1)/(exp(1)-1);
    for i=1:SearchAgents/2
        %% Phase1: Hunting and attacking strategy on iguana (Exploration Phase)
        iguana=Xbest;
        I=round(1+rand(1,1));

        X_P1(i,:)=w*X(i,:)+rand(1,1) .* (iguana-I.*X(i,:)); 
        X_P1(i,:) = max(X_P1(i,:),lowerbound);X_P1(i,:) = min(X_P1(i,:),upperbound);

        % update position
        L=X_P1(i,:);
        F_P1(i)=fitness(L);
        if(F_P1(i)<fit(i))
            X(i,:) = X_P1(i,:);
            fit(i) = F_P1(i);
        end

    end
    %%

    for i=1+SearchAgents/2 :SearchAgents

        iguana=lowerbound+rand(1,dimension).*(upperbound-lowerbound); %Eq(5)
        L=iguana;
        F_HL=fitness(L);
        I=round(1+rand(1,1));

        if fit(i)> F_HL
            X_P1(i,:)=w*X(i,:)+rand(1,1).*(iguana-I.*X(i,:)); %
        else
            X_P1(i,:)=w*X(i,:)+rand(1,1) .* (X(i,:)-iguana); % 
        end
        X_P1(i,:) = max(X_P1(i,:),lowerbound);X_P1(i,:) = min(X_P1(i,:),upperbound);

        % update position
        L=X_P1(i,:);
        F_P1(i)=fitness(L);
        if(F_P1(i)<fit(i))
            X(i,:) = X_P1(i,:);
            fit(i) = F_P1(i);
        end
    end
    %% END Phase1: Hunting and attacking strategy on iguana (Exploration Phase)

    %% Adaptive T-distribution variation strategy for process of escaping from predators%%
    freen = exp(4.*(t/Max_iterations).^2);  
    %% Phase2: The process of escaping from predators (Exploitation Phase)
    F_avg=mean(fit);
    for i=1:SearchAgents/2
        LO_LOCAL=lowerbound/t;% Eq(9)
        HI_LOCAL=upperbound/t;% Eq(10)
        if fit(i)> F_avg
            X_P2(i,:) =  Xbest+ trnd(freen)*Xbest;   
        else
            X_P2(i,:)=X(i,:)+(1-2*rand).* (LO_LOCAL+rand(1,1) .* (HI_LOCAL-LO_LOCAL)); % Eq. (8)
        end
        X_P2(i,:) = max(X_P2(i,:),LO_LOCAL);X_P2(i,:) = min(X_P2(i,:),HI_LOCAL);
        % update position based on Eq (11)
        L=X_P2(i,:);
        F_P2(i)=fitness(L);
        if(F_P2(i)<fit(i))
            X(i,:) = X_P2(i,:);
            fit(i) = F_P2(i);
        end
    end

    [fmax,idx]=max(fit);
    worse = X(idx,:);
    %% Alert mechanism for process of escaping from predators%%
    for i=1+SearchAgents/2 :SearchAgents
        if( fit(i)>fbest)
            X_P2(i,:)=Xbest+(randn(1,dimension)).*(abs(( X(i,:)-Xbest)));
        else
            X_P2(i,:) =X(i,:)+(2*rand(1)-1)*(abs(X(i,:)-worse))/ ( fit(i)-fmax+1e-50);
        end
        X_P2(i,:) = max(X_P2(i,:),lowerbound);X_P2(i,:) = min(X_P2(i,:),upperbound);
        % update position based on Eq (11)
        L=X_P2(i,:);
        F_P2(i)=fitness(L);
        if(F_P2(i)<fit(i))
            X(i,:) = X_P2(i,:);
            fit(i) = F_P2(i);
        end
    end

    %% END Phase2: The process of escaping from predators (Exploitation Phase)


    best_so_far(t)=fbest;
    average(t) = mean (fit);

end
Best_score=fbest;
Best_pos=Xbest;
TNTWCOA_curve=best_so_far;
end

