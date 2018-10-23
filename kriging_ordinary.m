function [f_out, std_out] = kriging_ordinary(sill, range, type, A, DS,DM)

Krigingvariableok=zeros(length(DS'),1);
Krigingstdok=zeros(length(DS'),1);

for index=1:length(DS')
    
    %At first, all variables need to be cleaned or reset because every 
    %iteration will modify them. 
    clear KDP k K Lambdask aux AUX DDP VV AUX2 excluderows npoints
    excluderows=0;
    aux=1;
    
    %This loop evaluates the distance between points and the "if"
    %segregates the ones, which are smaller than the chosen range
    for i=1:length(A')
%          a= distance(A(i,1),A(i,2),DS(index,1),DS(index,2),Earth); %Taking the Geodesic Distance
         a = sqrt((A(i,1)-DS(index,1))^2 + (A(i,2)-DS(index,2))^2); %In case we want to evaluate the Cartesian Distance
        if a<=range
        KDP(aux,1)=a;%KDP stands for Kriging Distance Points
        KDP(aux,2)=i;%This second column was created in order to know which points are being taken into account.
        aux=aux+1;
        end
    end
    
     %This step sorts the KDP matrix by the distance between points
     [AUX,I]=sort(KDP(:,1));
     KDP=KDP(I,:);
     
     %Using 16 points or KDP size, the minimun between them
     npoints=min(size(KDP,1),16);
     
     %The number '16' was chosen in order to get a considerable number of
     %points, but can be changed because it is somehow randomic.

    %Measuring distance between pairs of data points in order to exclude
    %the ones with distance greater than the range. Imagine that you have
    %3 aligned points, which we will call A,B and C. The point in the 
    %middle (B) is the one we want to apply the kriging process. And the
    %distances between AB and BC are smaller than the range. But, if the
    %distance between AC is greater than the range, wether A or C must be
    %deleted. In our case we are deleting the farest one. It is worth
    %highlighting that after excluding these problematic points it is 
    %necessary to check if are there another points that checks all boxes.
    %In case yes, these points must be taken into account.

    DDP=zeros(npoints);
    for i=1:npoints
        for j=1:npoints
        DDP(i,j)= DM(KDP(i,2),KDP(j,2)); %DDP stands for Distance of Data Points
        end
    end
    
    %This "if" verifies if is there any distance greater than the range
    %between the data points.
    if max(max(DDP))>range
        excluderows=1;
    else
    end
    
    %If any point must be excluded, then the program enters in the following
    %if.
    if excluderows==1
        clear AUX;
        a=1;

        %Note that before we were not interested in which point we shlould
        %exclude. We just decided that some point must be excluded.
        %This time we are looking for points that must be deleted. 
        
        while max(max(DDP))>range            
        %This 'while' command was used to garantee that at the end no
        %points with distance greater then the range will be taken.
        
            for i=1:npoints
                for j=1:npoints
                    if DDP(i,j)>range
                        AUX(1,1)=i;
                        AUX(1,2)=j;
                        break;
                    end
                end
            end

             %This 'if' excludes the farest point
             if KDP(AUX(1,1),1)>KDP(AUX(1,2),1)
                 AUX2=AUX(1,1);
             elseif KDP(AUX(1,1),1)<KDP(AUX(1,2),1)
                 AUX2=AUX(1,2);
             elseif KDP(AUX(1,1),1)==KDP(AUX(1,2),1)
                AUX2=AUX(1,2);
             end            
            
             %Removing Lines which contains undesired points
             KDP=removerows(KDP, AUX2);



             npoints=min(size(KDP,1),16);

             clear DDP;
             DDP=zeros(npoints);

            for i=1:npoints
                for j=1:npoints
                %DDP -> Distance of data points
                DDP(i,j)= DM(KDP(i,2),KDP(j,2));
                end
            end    
        end
    end
    
    %This next step calculates the semivariance using the chosen
    %semivariogram, which in this case is the Spherical Model. 

    k=zeros(npoints+1,2);%This matrix stores the semivariance as well the respective point.
    
    %It is worth saying that the semivariance measured in this loop is the
    %one between the point we are applying the kriging and the data points
    
    if strcmpi(type, 'exponential') ==1 
        for i=1:npoints     
            k(i,1)=sill*(1-exp((-3*KDP(i,1))/range));
            k(i,2)=KDP(i,2);
        end
    elseif strcmpi(type, 'spherical') ==1    
        for i=1:npoints    
        k(i,1)=sill*(1-1.5*(KDP(i,1)/range) + 0.5*(KDP(i,1)/range)^3);
        k(i,2)=KDP(i,2);
        end 
    elseif strcmpi(type, 'gaussian') ==1  
        for i=1:npoints     
            k(i,1)= sill*(1-exp((-3*power(KDP(i,1),2))/(range^2)));
            k(i,2)=KDP(i,2);
        end        
    else
        error;
    end
    
    %In this next loop we will calculate the distance of the data points
    %as well the semivariance between each other.
    
    K=zeros(npoints+1);  
    if strcmpi(type, 'exponential') ==1 
        for i=1:npoints
            for j=1:npoints
                DDP(i,j)= DM(k(i,2),k(j,2));
                K(i,j)= (1-exp((-3*DDP(i,j))/range));
            end
        end
    elseif strcmpi(type, 'spherical') ==1 
        for i=1:npoints
            for j=1:npoints
                DDP(i,j)= DM(k(i,2),k(j,2));
                K(i,j)=sill*(1-1.5*(DDP(i,j)/range) + 0.5*(DDP(i,j)/range)^3);
            end
        end
    elseif strcmpi(type, 'gaussian') ==1 
        for i=1:npoints
            for j=1:npoints
                DDP(i,j)= DM(k(i,2),k(j,2));
                K(i,j)=sill*(1-exp((-3*power(DDP(i,j),2))/(range^2)));
            end
        end
    else
        error;
    end
   
    %Until here all steps for Simple and Ordinary Kriging were exactly the
    %same. From now on it changes a little. As said before, the difference
    %between both krigings lies on the fact that the ordinary one assumes
    %that the mean is constant only in the local neighborhood. 
    %In matrix terms this difference reflects in the fact that this kriging
    %systemm is an augmented version of the simples kriging matrix. 
    %View -> http://people.ku.edu/~gbohling/cpe940 -> Kriging PDF (Pg. 15)
    %These next lines are responsible for writing the 'augmented matrix'
    
    for i=1:(npoints+1)
        if i<=npoints
            K(i,(npoints+1))= 1;
        else
        K(i,(npoints+1))=0;
        end
    end
    
    for i=1:npoints
            K((npoints+1), i)=1;
    end
    
    k((npoints+1),1)=1;
    
    Lambdaok= inv(K)*k(:,1); %Lambdaok stands for the lambda ordinary kriging vector
    %See -> http://people.ku.edu/~gbohling/cpe940 -> Kriging PDF

    for i=1:npoints
        VV(1,i)= A(k(i,2),3); %VV stands for Variable values
        %Once the mean (which is local) is included in the deduction of 
        %this method we do not have to subtract it. 
    end
    
    Lambdaok_ordinary=removerows(Lambdaok,(npoints+1));
    Krigingvariableok(index,1)= transpose(Lambdaok_ordinary)*VV';
    
    k=removerows(k, (npoints+1));
    aux= transpose(Lambdaok_ordinary)*k(:,1);
     Krigingstdok(index,1)=sqrt(abs(sill-aux-Lambdaok((npoints+1),1)));
          
end

f_out = Krigingvariableok;
std_out = Krigingstdok;
end