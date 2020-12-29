%RVEA Main File
function MAIN(Problem,M,Run)
%clc;
format compact;tic;%formate compact????????????

%basic settings
[Generations,N,p1,p2] = P_settings('RVEA',Problem,M);
Evaluations = Generations*N;
alpha = 2.0;
FE = 0;

%reference vector initialization
[N,Vs] = F_weight(p1,p2,M);
Vs(Vs==0) = 0.000000001;
for i = 1:N
    Vs(i,:) = Vs(i,:)./norm(Vs(i,:));
end;
V = Vs;
Generations = floor(Evaluations/N);

% %calculat neighboring angle for angle normalization
% cosineVV = V*V';
% [scosineVV, neighbor] = sort(cosineVV, 2, 'descend');%sort(A,2) sorts the elements of each row.
% acosVV = acos(scosineVV(:,2));%Inverse cosine in radians
% refV = (acosVV);

%population initialization
[Population,Boundary,Coding] = P_objective('init',Problem,M,N);
FunctionValue = P_objective('value',Problem,M,Population);

for Gene = 0 : Generations - 1
    %random mating and reproduction
    [MatingPool] = F_mating(Population);
    Offspring = P_generator(MatingPool,Boundary,Coding,N);  
    FE = FE + size(Offspring, 1);
    Population = [Population; Offspring];
    FunctionValue = [FunctionValue; P_objective('value',Problem,M,Offspring);];
    
    %APD based selection
    theta0 =  (Gene/(Generations))^alpha*(M);
    [Selection] = F_select(FunctionValue,V, theta0);
    Population = Population(Selection,:);
    FunctionValue = FunctionValue(Selection,:);

    %reference vector adaption
    if(mod(Gene, ceil(Generations*0.1)) == 0)
        %update the reference vectors
        Zmin = min(FunctionValue,[],1);	
        Zmax = max(FunctionValue,[],1);	
        V = Vs;
        V = V.*repmat((Zmax - Zmin)*1.0,N,1);
        for i = 1:N
            V(i,:) = V(i,:)./norm(V(i,:));
        end;
        %update the neighborning angle value for angle normalization
        cosineVV = V*V';
        [scosineVV, neighbor] = sort(cosineVV, 2, 'descend');
        acosVV = acos(scosineVV(:,2));
        refV = (acosVV); 
    end;

    clc; fprintf('Progress %4s%%\n',num2str(roundn(Gene/Generations*100,-1)));

end;
P_output(Population,toc,'RVEA',Problem,M,Run);
end


