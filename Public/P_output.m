function P_output (Population,time,Algorithm,Problem,M,Run)


FunctionValue = P_objective1('value',Problem,M,Population);
if(strcmp(Algorithm, 'cRVEA'))
    FunctionValue = FunctionValue(:,1:end - 1);
end;
TrueValue = P_objective1('true',Problem,M,1000);

NonDominated  = P_sort(FunctionValue,'first')==1;
Population    = Population(NonDominated,:);
FunctionValue1 = FunctionValue(NonDominated,:);
[HVvalue,PopObj] = HV(FunctionValue1,TrueValue);
IGDvalue = IGD(FunctionValue1,TrueValue);
PF=size(FunctionValue1,1);

if(M == 2)
    Plot2D(TrueValue, FunctionValue, 'ro');
end;

if(M == 3)
    Plot3D_PFboundry(TrueValue, FunctionValue, 'ro');
end;

eval(['save Data/',Algorithm,'/',Algorithm,'_',Problem,'_',num2str(M),'_',num2str(Run),' Population FunctionValue time HVvalue IGDvalue',])
end


