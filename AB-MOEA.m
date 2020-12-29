clc;
addpath RVEA;
addpath Public;
Problem='DTLZ6';%????
Algorithm = 'RVEA';
Pop_num=100;
objective=2;
Evaluations=199;
RunNum =1;
alpha=2;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%?????%%%%%%%%%%%%
tic
for Run=1:RunNum
    %% Generate random population
    M=objective;
    [Generations,D,p1,p2] = P_settings('RVEA',Problem,M);
    upboundry=ones(1,D);
    lpboundry=zeros(1,D);
    Boundary=[upboundry;lpboundry];
    N    = 11*D-1;
    L=11*D-1+25;
    THETA1 = 5.*ones(M,D);
    Model1 = cell(1,M);
    THETA2 = 5.*ones(M,D);
    Model2 = cell(1,M);
    Population = lhsamp(N,D);
    Population1=Population;
    Population2=Population;
    FunctionValue = P_objective1('value',Problem,M,Population);
    FunctionValue1=FunctionValue;
    FunctionValue2=FunctionValue;
    FE = size(Population,1);
    wmax=20;
    %% Optimization
    while FE<=Evaluations
        %%train GP
        % Delete duplicated solutions
        [~,index1]  = unique(Population1,'rows');
        PopDec1 = Population1(index1,:);
        PopObj1 = FunctionValue1(index1,:);
        Nondominate1 = P_sort(PopObj1,'first')==1;
        datalong1=size(PopDec1,1);        
        
        [~,index2]  = unique(Population2,'rows');
        PopDec2 = Population2(index2,:);
        PopObj2 = FunctionValue2(index2,:);
        Nondominate2 = P_sort(PopObj2,'first')==1;
        datalong2=size(PopDec2,1);
        %limit data
        if datalong1 <=L
            disp(sprintf('No training data decrease'));
            PopDec1=PopDec1;PopObj1=PopObj1 ;
        else
            paixu1 = non_domination_sort_mod([PopDec1, PopObj1],M);
            data1=paixu1(1:floor(L/2), 1:D+M);
            paixu1=paixu1(floor(L/2)+1:end, 1:D+M);
            index1=randperm(size(paixu1,1));
            data1=[data1;paixu1(index1(1:L-floor(L/2)), 1:D+M)];
            PopDec1=data1(:, 1:D);PopObj1=data1(:, D+1:D+M);
        end
         if datalong2 <=L
            disp(sprintf('No training data decrease'));
            PopDec2=PopDec2;PopObj2=PopObj2 ;
        else
            paixu2 = non_domination_sort_mod([PopDec2, PopObj2],M);
            data2=paixu2(1:floor(L/2), 1:D+M);
            paixu2=paixu2(floor(L/2)+1:end, 1:D+M);
            index2=randperm(size(paixu2,1));
            data2=[data2;paixu2(index2(1:L-floor(L/2)), 1:D+M)];
            PopDec2=data2(:, 1:D);PopObj2=data2(:, D+1:D+M);
        end
        for m=1:M
            % The parameter 'regpoly1' refers to one-order polynomial
            % function, and 'regpoly0' refers to constant function. The
            % former function has better fitting performance but lower
            % efficiency than the latter one
            dmodel1     = dacefit(PopDec1,PopObj1(:,m),'regpoly1','corrgauss',THETA1(m,:),1e-5.*ones(1,D),100.*ones(1,D));
            Model1{m}   = dmodel1;
            THETA1(m,:) = dmodel1.theta;
        end
         for m=1:M
            % The parameter 'regpoly1' refers to one-order polynomial
            % function, and 'regpoly0' refers to constant function. The
            % former function has better fitting performance but lower
            % efficiency than the latter one
            dmodel2     = dacefit(PopDec,PopObj(:,m),'regpoly1','corrgauss',THETA2(m,:),1e-5.*ones(1,D),100.*ones(1,D));
            Model2{m}   = dmodel2;
            THETA2(m,:) = dmodel2.theta;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%RVEA optimization
        w=1;
        [V0,Pop_num] = UniformPoint(Pop_num,M);
        V     = V0;
        V1    = V0;
        while w < wmax
            %             OffDec = GA(PopDec);
            [MatingPool] = F_mating(PopDec1);
            OffDec1 = P_generator(MatingPool,Boundary,'Real',N);
            PopDec1 = [PopDec1;OffDec1];
            [N,~]  = size(PopDec1);
            Popmean1 = zeros(N,M);
            MSE1    = zeros(N,M);
            for i = 1: N
                for j = 1 : M
                    [PopObj1(i,j),~,MSE1(i,j)] = predictor(PopDec1(i,:),Model1{j});
                end
            end
            Selection = FSelection(PopObj1,V,(w/wmax)^alpha);
            PopDec1 = PopDec1(Selection,:);
            PopObj1 = PopObj1(Selection,:);
            MSE1=MSE1(Selection,:);
            if(mod(w, ceil(wmax*0.1)) == 0)
                V = V0.*repmat(max(FunctionValue1,[],1)-min(FunctionValue1,[],1),size(V0,1),1);
            end;
            w=w+1;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%RVEA optimization
        w=1;
        [V0,Pop_num] = UniformPoint(Pop_num,M);
        V     = V0;
        V1    = V0;
        while w < wmax
            %             OffDec = GA(PopDec);
            [MatingPool] = F_mating(PopDec2);
            OffDec2 = P_generator(MatingPool,Boundary,'Real',N);
            PopDec2 = [PopDec2;OffDec2];
            [N,~]  = size(PopDec2);
            Popmean2 = zeros(N,M);
            MSE2    = zeros(N,M);
            for i = 1: N
                for j = 1 : M
                    [PopObj2(i,j),~,MSE2(i,j)] = predictor(PopDec2(i,:),Model2{j});
                end
            end
            Selection = FSelection(PopObj2,V,(w/wmax)^alpha);
            PopDec2 = PopDec2(Selection,:);
            PopObj2 = PopObj2(Selection,:);
            MSE2=MSE2(Selection,:);
            if(mod(w, ceil(wmax*0.1)) == 0)
                V = V0.*repmat(max(FunctionValue2,[],1)-min(FunctionValue2,[],1),size(V0,1),1);
            end;
            w=w+1;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%acquisiion function
        a=-0.5*cos(FE*pi/300)+0.5;

        if a>0.5
            F=2;
        else
            F=1;
        end
        b=0.5*cos(FE*pi/300)+0.5;
        
        [MMSE1,~]=max(MSE1,[],1);
        [MPopObj1,~]=max(PopObj1,[],1);
        [MMSE2,~]=max(MSE2,[],1);
        [MPopObj2,~]=max(PopObj2,[],1);
        
        fit1=PopObj1./MPopObj1*b+MSE1./MMSE1*a;
       
        %%%%%%%%%%%select by the reference vectors
        if F==2
            Selection1 = FSelection(fit1,V1,(FE/300)^alpha);
            PopNew1 = PopDec1(Selection1,:);
            if size(PopNew1,1)<5
                PopNew3=PopNew1(:,1:D);
            elseif size(PopNew1,1)>=5
                PopNew3=PopNew1(randperm(size(PopNew1,1),5),1:D);
            end
            Population1 = [Population1;PopNew3];
            
        end
        
        
        if F==1
             Selection = FSelection1(fit,V1,(FE/300)^alpha);
             PopNew1 = PopDec(Selection,:);
             if size(PopNew1,1)<5
                PopNew=PopNew1(:,1:D);
				elseif size(PopNew1,1)>=5
                 PopNew=PopNew1(randperm(size(PopNew1,1),5),1:D);
             end
         end
         if(mod(FE, ceil(300*0.1)) == 0)
            V1 = V0.*repmat(max(FunctionValue,[],1)-min(FunctionValue,[],1),size(V0,1),1);
        end;
        
        
     Population1 = [Population1;PopNew3];
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         Population2 = [Population2;PopNew4];
        FunctionValue1=[P_objective1('value',Problem,M,Population1)]
        FE = FE+5;
%         P_output(Population,toc,'RVEA_SM',Problem,M,Run);
    end 
        P_output1(Population1,'RVEA_SM',Problem,M);
        
    Run=Run+1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%when EI is a multi-objective value
%         candidate=[PopDec PopObj fit];
%         %%%%%%%%%%%%%%%%%%%select the PF by decision vector(bad)%%%%%%
%         candidateSort=non_domination_sort_mod(candidate,M);
%          Nondominate = find(candidateSort(:,D+2*M+1)==1);
%          if size(Nondominate,1)<=10
%              PopNew=candidate(randperm(size(Nondominate,1),min(size(Nondominate,1),5)),1:D);
%          elseif size(Nondominate,1)>15
%              PopNew= candidate(Nondominate,:);
%              cindex  = kmeans(PopNew(:,1:D),5);
%              Q = [];
%              for i = 1 : 5
%                  index = find(cindex == i);
%                  tempFit=PopNew(index,D+4:D+6);
%                  [Y,I]=min(mean(tempFit,2));
%                  [~,best] = min(Y);
%                  Q = [Q,index(best)];
%              end
%              PopNew = PopNew(Q,1:D);
%          end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%select PF by fit(better than decision vector)
%         candidateSort=non_domination_sort_mod(candidate,M);
%          Nondominate = find(candidateSort(:,D+2*M+1)==1);
%          if size(Nondominate,1)<=5
%              PopNew=candidate(randperm(size(Nondominate,1),min(size(Nondominate,1),5)),1:D);
%          elseif size(Nondominate,1)>5
%              PopNew= candidate(Nondominate,:);
%              cindex  = kmeans(PopNew(:,D+4:D+6),5);
%              Q = [];
%              for i = 1 : 5
%                  index = find(cindex == i);
%                  tempFit=PopNew(index,D+4:D+6);
%                  [Y,I]=min(mean(tempFit,2));
%                  [~,best] = min(Y);
%                  Q = [Q,index(best)];
%              end
%              PopNew = PopNew(Q,1:D);
%          end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%when EI is a single value(bad)
%         [~,numindex]=min(fit);
%         PopNew=PopDec(numindex,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%select randomly
%          if size(Nondominate,1)<5
%              PopNew=candidate(:,1:D);
%          elseif size(Nondominate,1)>=5
%             PopNew=candidate(randperm(size(Nondominate,1),5),1:D);
%          end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%the best currently
%         %%%%%%%%%%%%%%%%%%%select the PF%%%%%%
%         candidateSort=non_domination_sort_mod(candidate,M);
%          Nondominate = find(candidateSort(:,D+2*M+1)==1);
%          if size(Nondominate,1)<5
%              PopNew=candidate(:,1:D);
%          elseif size(Nondominate,1)>=5
%             PopNew=candidate(randperm(size(Nondominate,1),5),1:D);
%          end
%