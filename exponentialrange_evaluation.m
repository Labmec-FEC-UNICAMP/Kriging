function [Exponentialmodel_range,r2exponential] = exponentialrange_evaluation(DM,LAG, sill, SVV, Weight )

%This step creates the DSX vector
DSX=linspace(0, max(max(DM(:))), ((max(max(DM(:))))+1000));

%This loop replaces several times the LAG variable of the Gaussian Model 
%with a  value (from DSX).  After each replacement the loop measures how 
%good this lag was. 
%This measurement is done twice. The first one is done using a
%weighted nonlinear regression. The weightning is done by the number of
%data pairs for each lag. These values were stored by the 'weight' vector.
%The other method used to determine wether the lag is good or not is the
%'r-squared' method. In all tested cases the output of each method was the
%same. 

minerror=false;
r2exponential=false;
Exponentialmodel=zeros(length(DSX),1);

for i=1:length(DSX)

    aux=0;
    aux2=0;
    aux3=0;
    for j=1:length(LAG)
    Exponentialmodel(i,1)=sill*(1-exp((-3*LAG(1,j))/(DSX(1,i) + eps)));
    aux=((Exponentialmodel(i,1)-SVV(j,1))^2)*(1/Weight(j,1));
    aux2=aux2 + aux;
    clear aux
    aux=((SVV(j,1)-mean(SVV))^2);
    aux3=aux3 + aux;
    end
    if aux2<minerror ||minerror == false
        minerror=aux2;
        Exponentialmodel_range=DSX(1,i);
    end
    rsquaredexponential=1-(aux2/aux3);
    if rsquaredexponential>r2exponential && rsquaredexponential>=0 && rsquaredexponential<=1 || r2exponential ==false
        r2exponential=rsquaredexponential;
        Exponentialmodel_range2=DSX(1,i);
    end
end


Exponentialmodel_range = Exponentialmodel_range;
r2exponential = r2exponential;

end