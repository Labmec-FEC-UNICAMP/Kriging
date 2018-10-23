function sill = sill_evaluation(A, LAGM)

k=1;

for i=1:length(LAGM');

        if LAGM(i,1)==0
        TAIL(k,1)= A(LAGM(i,2),3);
        HEAD(k,1)= A(LAGM(i,3),3);
        k=k+1;
        end
        
end

%Calculating the Covariance
mtail=mean(TAIL);
mhead=mean(HEAD);
aux=0;

for i=1:k-1
    aux =aux+ TAIL(i,1)*HEAD(i,1);
end

aux=aux/(k-1);
sill=aux-mtail*mhead;

sill = sill;
end