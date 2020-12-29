function P_output1 (Population,Algorithm,Problem,M)
[r,l]=size(Population);
FunctionValue=P_objective1('value',Problem,M,Population);
FunctionValue1 = FunctionValue(r-9:r,:);

 NonDominated  = P_sort(FunctionValue(1:r-10,:),'first')==1;
 Population    = Population(NonDominated,:);
 FunctionValue2 = FunctionValue(NonDominated,:);
 
TrueValue = P_objective1('true',Problem,M,1000);

PF=size(FunctionValue2,1);

% if(M == 2)
%     Plot2D(TrueValue, FunctionValue, 'ro');
% end;

if(M == 3)
    Plot3D1_PFboundry(TrueValue, FunctionValue1,FunctionValue2, 'ro');
end;

% eval(['save Data/',Algorithm,'/',Algorithm,'_',Problem,'_',num2str(M),'_',num2str(Run),' Population FunctionValue time HVvalue IGDvalue',])
end


