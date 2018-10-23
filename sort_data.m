function [] = sort_data(A)

%Sorting A matrix by characteristic for interpolation
%It is also important to set that Z is the characteristic we are looking to
%interpolate
[~,I]=sort(A(:,3));
B=A(I,:);
clear I

%Splitting B matrix in 5 groups in order to plot the Histogram and also
%plotting them in space to check if there is some significant trend

%This step evaluate each group limits
ZV=linspace(min(B(:,3)),max(B(:,3)), 6); %ZV stands for Z values

%This loop splits all the groups, but it also returns some zeros
%Therefore will be necessary some lines of code to discard all these
%unused numbers

XV=NaN([length(B') length(ZV)]);
YV=NaN([length(B') length(ZV)]);
VV=NaN([length(B') length(ZV)]);



for j=1:length(ZV)
    aux=0;
    for i=1:length(B')        
        if j<length(ZV)&& A(i,3)>= ZV(1,j) && A(i,3)<= ZV(1,j+1)
            aux=aux+1;
            XV(aux,j) = A(i,2);
            YV(aux,j) = A(i,1);
            VV(aux,j) = A(i,3);
        end   
    end  
end

%Then it will split the whole matrix into groups and also clear all lines
%which contians only zeros
XV1=XV(:,1);
XV1(isnan(XV1)) = [];
XV2=XV(:,2);
XV2(isnan(XV2)) = [];
XV3=XV(:,3);
XV3(isnan(XV3)) = [];
XV4=XV(:,4);
XV4(isnan(XV4)) = [];
XV5=XV(:,5);
XV5(isnan(XV5)) = [];

YV1=YV(:,1);
YV1(isnan(YV1)) = [];
YV2=YV(:,2);
YV2(isnan(YV2)) = [];
YV3=YV(:,3);
YV3(isnan(YV3)) = [];
YV4=YV(:,4);
YV4(isnan(YV4)) = [];
YV5=YV(:,5);
YV5(isnan(YV5)) = [];

%Note that the groups were split by wind speed, BUT the speed itself is
%not necessary to the plot. Only its longitude and latitude are enough for
%this plot
figure
P1 = scatter(XV1,YV1,'y','filled');
hold on;
P2 = scatter(XV2,YV2,'m','filled');
hold on;
P3 = scatter(XV3,YV3,'c','filled');
hold on;
P4 = scatter(XV4,YV4,'r','filled');
hold on;
P5 = scatter(XV5,YV5,'g','filled');


%Setting Graph parameters
% title('Wind Speed (m/s)')
% xlabel('Longitude (Decimal Systemm)')
% ylabel('Latitude (Decimal Systemm)')

%In order to create an organized plot, with legends the borders of each
%group will be converted into strings
LG1=strcat(num2str(ZV(1,1),4), '-' , num2str(ZV(1,2),4));
LG2=strcat(num2str(ZV(1,2),4), '-' , num2str(ZV(1,3),4));
LG3=strcat(num2str(ZV(1,3),4), '-' , num2str(ZV(1,4),4));
LG4=strcat(num2str(ZV(1,4),4), '-' , num2str(ZV(1,5),4));
LG5=strcat(num2str(ZV(1,5),4), '-' , num2str(ZV(1,6),4));

%Setting legend parameters
legend([P1 P2 P3 P4 P5],LG1,LG2,LG3,LG4,LG5);
legend('Location', 'southeast')
legend('Boxoff')

end