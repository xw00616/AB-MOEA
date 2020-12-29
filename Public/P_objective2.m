%This file includes the original and scaled DTLZ1 to DTLZ4

function [Output,Boundary,Coding] = P_objective2(Operation,Problem,M,Input)
k = find(~isstrprop(Problem,'digit'),1,'last');
switch Problem(1:k)
    case 'UF'
        [Output,Boundary,Coding] = P_UF(Operation,Problem,M,Input);
    case 'SDTLZ'
        [Output,Boundary,Coding] = P_SDTLZ(Operation,Problem,M,Input);
    otherwise
        error([Problem,'Not Exist']);
end
end

function [Output,Boundary,Coding] = P_UF(Operation,Problem,M,Input)
persistent K;
Boundary = NaN; Coding = NaN;
switch Operation
    %Population Initialization
    case 'init'
        D = 30;
        N=Input;
        switch Problem
            case  {'UF1','UF2','UF5','UF6','UF7'}
                lower    = [0,zeros(1,D-1)-1];
                upper    = ones(1,D);
                Population= unifrnd(repmat(lower,N,1),repmat(upper,N,1));
%                 Population= repmat(MinValue,Input,1) + repmat(MaxValue-MinValue,Input,1).*lhsdesign(Input,D,'criterion','maximin','iterations',1000);
                Output   = Population;
                Boundary = [lower;upper];
                Coding   = 'Real';
            case 'UF3'
                lower    = zeros(1,D);
                upper    = ones(1,D);
                 Population= unifrnd(repmat(lower,N,1),repmat(upper,N,1));
                Output   = Population;
                Boundary = [lower;upper];
                Coding   = 'Real';
            case 'UF4'
                lower    = [0,zeros(1,D-1)-2];
                upper    = [1,zeros(1,D-1)+2];
                 Population= unifrnd(repmat(lower,N,1),repmat(upper,N,1));
                Output   = Population;
                Boundary = [lower;upper];
                Coding   = 'Real';
                %Objective Function Evaluation
        end
            case 'value'
                X   = Input;
                PopObj = zeros(size(X,1),M);
                switch Problem
                    case 'UF1'
                        D  = size(X,2);
                        J1 = 3 : 2 : D;
                        J2 = 2 : 2 : D;
                        Y  = X - sin(6*pi*repmat(X(:,1),1,D)+repmat(1:D,size(X,1),1)*pi/D);
                        PopObj(:,1) = X(:,1)         + 2*mean(Y(:,J1).^2,2);
                        PopObj(:,2) = 1-sqrt(X(:,1)) + 2*mean(Y(:,J2).^2,2);
                        
                    case 'UF2'
                        D  = size(X,2);
                        J1 = 3 : 2 : D;
                        J2 = 2 : 2 : D;
                        Y       = zeros(size(X));
                        X1      = repmat(X(:,1),1,length(J1));
                        Y(:,J1) = X(:,J1)-(0.3*X1.^2.*cos(24*pi*X1+4*repmat(J1,size(X,1),1)*pi/D)+0.6*X1).*cos(6*pi*X1+repmat(J1,size(X,1),1)*pi/D);
                        X1      = repmat(X(:,1),1,length(J2));
                        Y(:,J2) = X(:,J2)-(0.3*X1.^2.*cos(24*pi*X1+4*repmat(J2,size(X,1),1)*pi/D)+0.6*X1).*sin(6*pi*X1+repmat(J2,size(X,1),1)*pi/D);
                        PopObj(:,1) = X(:,1)         + 2*mean(Y(:,J1).^2,2);
                        PopObj(:,2) = 1-sqrt(X(:,1)) + 2*mean(Y(:,J2).^2,2);
                    case 'UF3'
                        D  = size(X,2);
                        J1 = 3 : 2 : D;
                        J2 = 2 : 2 : D;
                        Y  = X - repmat(X(:,1),1,D).^(0.5*(1+3*(repmat(1:D,size(X,1),1)-2)/(D-2)));
                        PopObj(:,1) = X(:,1)         + 2/length(J1)*(4*sum(Y(:,J1).^2,2)-2*prod(cos(20*Y(:,J1)*pi./sqrt(repmat(J1,size(X,1),1))),2)+2);
                        PopObj(:,2) = 1-sqrt(X(:,1)) + 2/length(J2)*(4*sum(Y(:,J2).^2,2)-2*prod(cos(20*Y(:,J2)*pi./sqrt(repmat(J2,size(X,1),1))),2)+2);
                    case 'UF4'
                        D  = size(X,2);
                        J1 = 3 : 2 : D;
                        J2 = 2 : 2 : D;
                        Y  = X - sin(6*pi*repmat(X(:,1),1,D)+repmat(1:D,size(X,1),1)*pi/D);
                        hY = abs(Y)./(1+exp(2*abs(Y)));
                        PopObj(:,1) = X(:,1)      + 2*mean(hY(:,J1),2);
                        PopObj(:,2) = 1-X(:,1).^2 + 2*mean(hY(:,J2),2);
                    case 'UF5'
                        D  = size(X,2);
                        J1 = 3 : 2 : D;
                        J2 = 2 : 2 : D;
                        Y  = X - sin(6*pi*repmat(X(:,1),1,D)+repmat(1:D,size(X,1),1)*pi/D);
                        hY = 2*Y.^2 - cos(4*pi*Y) + 1;
                        PopObj(:,1) = X(:,1)   + (1/20+0.1)*abs(sin(20*pi*X(:,1)))+2*mean(hY(:,J1),2);
                        PopObj(:,2) = 1-X(:,1) + (1/20+0.1)*abs(sin(20*pi*X(:,1)))+2*mean(hY(:,J2),2);
                    case 'UF6'
                        D  = size(X,2);
                        J1 = 3 : 2 : D;
                        J2 = 2 : 2 : D;
                        Y  = X - sin(6*pi*repmat(X(:,1),1,D)+repmat(1:D,size(X,1),1)*pi/D);
                        PopObj(:,1) = X(:,1)   + max(0,2*(1/4+0.1)*sin(4*pi*X(:,1)))+2/length(J1)*(4*sum(Y(:,J1).^2,2)-2*prod(cos(20*Y(:,J1)*pi./sqrt(repmat(J1,size(X,1),1))),2)+2);
                        PopObj(:,2) = 1-X(:,1) + max(0,2*(1/4+0.1)*sin(4*pi*X(:,1)))+2/length(J2)*(4*sum(Y(:,J2).^2,2)-2*prod(cos(20*Y(:,J2)*pi./sqrt(repmat(J2,size(X,1),1))),2)+2);
                        
                    case 'UF7'
                        D  = size(X,2);
                        J1 = 3 : 2 : D;
                        J2 = 2 : 2 : D;
                        Y  = X - sin(6*pi*repmat(X(:,1),1,D)+repmat(1:D,size(X,1),1)*pi/D);
                        PopObj(:,1) = X(:,1).^0.2   + 2*mean(Y(:,J1).^2,2);
                        PopObj(:,2) = 1-X(:,1).^0.2 + 2*mean(Y(:,J2).^2,2);
                        
                end
                Output =  PopObj;
                %Sample True PFs
            case 'true'
                N=Input;
                switch Problem
                    
                    case {'UF1','UF2','UF3'}
                        P(:,1) = (0:1/(N-1):1)';
                        P(:,2) = 1 - P(:,1).^0.5;
                    case 'UF4'
                        P(:,1) = (0:1/(N-1):1)';
                        P(:,2) = 1 - P(:,1).^2;
                    case 'UF5'
                        P(:,1) = (0:1:20)'/20;
                        P(:,2) = 1 - P(:,1);
                    case 'UF6'
                        P(:,1) = (0:1/(N-1):1)';
                        P(:,2) = 1 - P(:,1);
                        P(P(:,1)>0 & P(:,1)<1/4 | P(:,1)>1/2 & P(:,1)<3/4,:) = [];
                    case 'UF7'
                        P(:,1) = (0:1/(N-1):1)';
                        P(:,2) = 1 - P(:,1);
                end
                Output = P;
 
end
end



function [Output,Boundary,Coding] = P_SDTLZ(Operation,Problem,M,Input)
persistent K;
persistent F;
Boundary = NaN; Coding = NaN;
switch Operation
    %Population Initialization
    case 'init'
        k = find(~isstrprop(Problem,'digit'),1,'last');
        K = [5 10 10 10 10 10 20];
        K = K(str2double(Problem(k+1:end)));
        F = [10 10 10 10 10 5 4 3 2 2];
        F = F(M);
        
        D = M+K-1;
        MaxValue   = ones(1,D);
        MinValue   = zeros(1,D);
        Population = rand(Input,D);
        Population = Population.*repmat(MaxValue,Input,1)+(1-Population).*repmat(MinValue,Input,1);
        
        Output   = Population;
        Boundary = [MaxValue;MinValue];
        Coding   = 'Real';
        %Objective Function Evaluation
    case 'value'
        Population    = Input;
        FunctionValue = zeros(size(Population,1),M);
        switch Problem
            case 'SDTLZ1'
                g = 100*(K+sum((Population(:,M:end)-0.5).^2-cos(20.*pi.*(Population(:,M:end)-0.5)),2));
                for i = 1 : M
                    FunctionValue(:,i) = 0.5.*prod(Population(:,1:M-i),2).*(1+g);
                    if i > 1
                        FunctionValue(:,i) = FunctionValue(:,i).*(1-Population(:,M-i+1));
                    end
                end
            case 'SDTLZ2'
                g = sum((Population(:,M:end)-0.5).^2,2);
                for i = 1 : M
                    FunctionValue(:,i) = (1+g).*prod(cos(0.5.*pi.*Population(:,1:M-i)),2);
                    if i > 1
                        FunctionValue(:,i) = FunctionValue(:,i).*sin(0.5.*pi.*Population(:,M-i+1));
                    end
                end
            case 'SDTLZ3'
                g = 100*(K+sum((Population(:,M:end)-0.5).^2-cos(20.*pi.*(Population(:,M:end)-0.5)),2));
                for i = 1 : M
                    FunctionValue(:,i) = (1+g).*prod(cos(0.5.*pi.*Population(:,1:M-i)),2);
                    if i > 1
                        FunctionValue(:,i) = FunctionValue(:,i).*sin(0.5.*pi.*Population(:,M-i+1));
                    end
                end
            case 'SDTLZ4'
                Population(:,1:M-1) = Population(:,1:M-1).^100;
                g = sum((Population(:,M:end)-0.5).^2,2);
                for i = 1 : M
                    FunctionValue(:,i) = (1+g).*prod(cos(0.5.*pi.*Population(:,1:M-i)),2);
                    if i > 1
                        FunctionValue(:,i) = FunctionValue(:,i).*sin(0.5.*pi.*Population(:,M-i+1));
                    end
                end
        end
        Output = FunctionValue.*repmat((F.^(0:M - 1)), [size(FunctionValue,1) 1]);
        %Sample True PFs
    case 'true'
        switch Problem
            case 'SDTLZ1'
                Population = T_uniform(Input,M);
                Population = Population/2;
            case {'SDTLZ2','SDTLZ3','SDTLZ4'}
                Population = T_uniform(Input,M);
                for i = 1 : size(Population,1)
                    Population(i,:) = Population(i,:)./norm(Population(i,:));
                end
        end
        Output = Population.*repmat((F.^(0:M - 1)), [size(Population,1) 1]);
end
end

