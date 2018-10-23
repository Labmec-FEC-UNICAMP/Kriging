function [Gaussianmodel_range, r2gaussian] = gaussianrange_evaluation(DM,LAG, sill, SVV, Weight )

%This step creates the DSX vector
DSX=linspace(0, max(max(DM(:))), ((max(max(DM(:))))+1));

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
r2gaussian=false;
Gaussianmodel=zeros(length(DSX),1);
%Minerror stands for minimun error and r2 to r-squared. At first, they were 
%set as 'false' to ensure that it gets in the 'if' function. But, later 
%they will become  a number.

for i=1:length(DSX)
    
    %These are auxiliar variables and must be zeroed every itaration. 
    aux=0;
    aux2=0;
    aux3=0;
    
    
    for j=1:length(LAG)
    Gaussianmodel(i,1)=sill*(1-exp((-3*power(LAG(1,j),2))/(DSX(1,i)^2 + eps)));
    aux=((Gaussianmodel(i,1)-SVV(j,1))^2)*(1/Weight(j,1));
    aux2=aux2 + aux;
    clear aux
    aux=((SVV(j,1)-mean(SVV))^2);
    aux3=aux3 + aux;
    end
    
    if aux2<minerror ||minerror == false
        minerror=aux2;
        Gaussianmodel_range=DSX(1,i);
    end
    
    rsquaredgaussian=1-(aux2/aux3);
    
    if rsquaredgaussian>r2gaussian && rsquaredgaussian>=0 && rsquaredgaussian<=1 || r2gaussian ==false
        r2gaussian=rsquaredgaussian;
        Gaussianmodel_range2=DSX(1,i);
    end
    
end

Gaussianmodel_range = Gaussianmodel_range;
r2gaussian = r2gaussian;

end
