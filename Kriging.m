clc;
clear all;
close all;
%% Input Data
%% Importing text file 

A = dlmread('TesteKriging_Parametrizado.txt');

%Data info
%Latitude   (Decimal)           Distance North or South of the Equator
%Longitude  (Decimal)           Distance East or West of the Greenwich 
%Wind Basic Speed (m/s)

%% Comparing the input data with normal distribution
%% Searching For Significant Trends

%Sorting A matrix by Z charachteristic
%A must be in the correct format. It means: 1st column -> x; 2nd column ->
%y; 3rd column z (the variable we want to predict)
sort_data(A);

%% Doing the Histogram
%In order to plot an histogram superimposed by a normal curve, an x vector
%will be discretized  and some other needed parameters will be
%calculated

x = min(A(:,3)):0.01:max(A(:,3));
meanVV = mean(A(:,3));                  %meanVV stands for mean of the Variable Vector
stdVV = std(A(:,3));                    %stdVV stands for standard deviaton of the Variable Vector 

%Doing the histogram
figure
h = histogram(A(:,3),'Normalization','pdf');

%It's necessary to set how many columns will have this histogram
%In this case was set 6 columns
h.BinWidth = (max(A(:,3)) - min(A(:,3)))/5;
hold on
plot(x,normpdf(x,meanVV,stdVV))

%Setting labels
xlabel('Z')
ylabel('Relative Frequency')

%% Normal Quantile Plot

%In Matlab is possible to find a function which does all the QQ-Plot
%needed calculations. Therefore, it only misses to plot.
figure
qqplot(A(:,3))

%% KS - Test

%Transform data into standard normal results
x = (A(:,3)-meanVV)/(stdVV);

%Null-hypothesis (is this data well represented by a normal distribution?)
%If h = 1, this indicates the rejection of the null hypothesis at the Alpha
%significance level.
%If h = 0, this indicates a failure to reject the null hypothesis at the
%Alpha significance level. The returned value of h = 0 indicates that
%kstest fails to reject the null hypothesis at the default 5% significance
%level


[h,p] = kstest(x);
%what is the p-value? p-value of the test, returned as a scalar value in
%the range [0,1]. Small values of p cast doubt on the validity of the null
%hypothesis.

%Create a CDF based on the empirical results (x)
[f,x_values] = ecdf(x);

figure
P1 = plot(x_values,f);
hold on;

%Plot a standard normal CDF (mean=0 and std=1)
P2 = plot(x_values,normcdf(x_values,0,1),'r-');

%Setting Graph parameters
title('Cumulative Distribution Function')
xlabel('Standard Normal Quantiles')
ylabel('Cumulative Probability')

%Setting legend parameters
legend([P1 P2], 'Empirical CDF','Normal CDF');
legend('Location', 'southeast')
legend('Boxoff')

 %% Important Checkpoint!!
 
 %Untill this point the program was doing tests in order to compare wether
 %the input data was similar enough to a normal distribution or not.Please
 %note that if the input data is not similar enough to a normal 
 %distribution the kriging process should not be executed
 
 %%
 %% Distance Measurement
 
%The kriging process is all about distances. Therefore, to simplify the
%future calculations a Distance Matrix (DM) will be evaluated. It will
%contains the distance between every input point

%It's also important to measure the right distance. In some cases the
%cartesian coordinates are enough, but in some cases it's necessary to
%evaluate the geodesic distance. The function below evaluete both distances
%types. You just need to set which distance type is better for you.

DM = distance_measurement(A, 'cartesian');

%%
%% Plotting Statistics versus Lag
%As said before, the greater is the lag, the smaller is the correlation
%between each point. In order to let it clear, this section intends to
%create a graph showing the behaviour of correlation, covariance and the 
%semivariance when varying the lag. It was possible to see it in the 
%crossplots. Now, the same things are going to be done, however, with 
%different lags. Note that the range is still the same.

%The first step is to reshape the DM (Distance Matrix) into a columm matrix
%The second and third columns were set to ease the understanding of each
%line. This vector will be called LAGM.

LAGM = zeros((0.5*((length(A'))^2 + length(A'))),3);
%LAGM stands for Lag Matrix

k = 1;
for i = 1:length(A')
    for j = 1:length(A')
        if i <= j
        LAGM(k,1) = DM(i,j);
        LAGM(k,2) = i;
        LAGM(k,3) = j;
        k = k+1;
        end
    end
end



%These are the LAGs 
spc = 1000; %space between lags
LAG = linspace(0, 0.5*(roundn(max(LAGM(:,1)), round(log10(max(LAGM(:,1)))))), ((0.5*(roundn(max(LAGM(:,1)), round(log10(max(LAGM(:,1)))))))/spc +1));
% This formula rounds the lags to distances like 10,20 or 100,200 or
% 1000,2000 and so on. 

%This net loop will evaluate the Covariance, Semivariance 
%and the Correlation for each lag in the vector 'LAG'. After, it will chart
%the calculated parameters versus its lags. The point of this loop is to
%observe the parameters varying versus the lag.
%The logical output is a decreasing covariance and correlation and a
%increasing semivariance.

Weight = zeros(length(LAG),1);
SVV = zeros(length(LAG),1);
COVV = zeros(length(LAG),1);
CORRV = zeros(length(LAG),1);

for i = 1:length(LAG)
   k = 1;
   clear TAIL;
   clear HEAD;
    for j = 1:length(LAGM');
        if LAGM(j,1) >= (LAG(1,i)-(max(LAGM(:,1))/36)) && (LAGM(j,1)<=LAG(1,i)+(max(LAGM(:,1))/36))
            TAIL(k,1) = A(LAGM(j,2),3);
            HEAD(k,1) = A(LAGM(j,3),3);
            k=k+1;
        end
    end
    
    %The 'Weight' vector will be further used. It stores how many points
    %each lag has
    Weight(i,1) = k-1;
     
    %Calculating the semivariance
    aux = 0;
    for a = 1:k-1
        aux = aux + (HEAD(a,1) - TAIL(a,1))^2;
    end
    SVV(i,1) = aux/(2*(k-1));
    
    %Calculating the covariance
    mtail = mean(TAIL);
    mhead = mean(HEAD);
    aux=0;
    for a = 1:k-1
        aux = aux+ TAIL(a,1)*HEAD(a,1);
    end
    aux = aux/(k-1);
    COVV(i,1) = aux-mtail*mhead;
    
    %Calculating the correlation
    aux = 0;
    for a = 1:k-1
        aux = aux + (TAIL(a,1)-mtail)^2;
    end
    sdtail = aux/(k-1);

    aux = 0;
    for a = 1:k-1
        aux = aux + (HEAD(a,1)-mtail)^2;
    end
    sdhead = aux/(k-1);

    CORRV(i,1) = COVV(i,1)/(sqrt(sdtail*sdhead));
end

%Setting Graph Parameters
figure
plot (LAG, SVV,'r--o');
hold on;
plot (LAG, COVV, 'b:+');

%Setting legend parameters
title('Semivariance and Covariance versus Mean Lag')
xlabel('Lag (km)')
legend('Semivariance', 'Covariance')
legend('Location', 'northwest')
legend('Boxoff')

figure
plot (LAG,CORRV, 'g-.*');
hold on;
title('Correlation versus Mean Lag')
xlabel('Lag (km)')
legend('Correlation')
legend('Boxoff')

%% Evaluating the Semivariogram parameters
%First of all it is important to understand what a semivariogram
%represents. Theoretically, the semivariogram is a curve which represents
%every semivariance versus every lag. So, the exact semivariogram could be
%obtained by calculating every semivariance for every lag. But, it would
%not result in a smooth curve and would take too much time. Beyond, the
%the kriging process demands a smooth semivariogram.
%Usually, to work around this problem, some semivariance x LAG are plotted
%in a graph and then an adjustment curve is set. This adjustment curve is
%called 'semivariogram'.
%In order to plot a semivariogram three basics parameters are needed. The
%first one is called 'Sill'. The sill is the value which the semivariogram
%flattens itself. The second one is the 'Range'. The range consists of the
%value which the semivariogram reaches the sill. The third one is the
%nugget. This point means the beginning of the semivariogram.
%It's worth remebering that the nugget should always be null. However,
%not always the calculated is equal to zero. In this program the nugget is
%not evaluated and is forced to be null. 

%% Sill
%The first parameter to be calculated is the 'Sill'. The sill corresponds
%to the covariance when the lag is exactly 0. In other words, it can be
%said that the sill is the covariance of all input points.
%Thus, to measure the sill is needed to determine the covariance, the same
%way it was done before. 
% clear TAIL;
% clear HEAD;

sill = sill_evaluation(A, LAGM);


%% Range
%As said before, the range corresponds to the lag distance which the
%semivariogram flattens itself. Theoretically, this is the turning
%point of the whole process. For any lag greater than the range, the
%correlation is essentially zero. At first, it seems a good idea to
%calculate the correlation for several lags and then pick the lag when
%correlation is null. However, this is not a good idea. There are some
%reasons for not doing that. The fisrt one is that there is not any
%difference between correlation 0.1 or 0 or -0.1, i.e., which one is null
%enough? Why is 0 better or worse than 0.1. In matter of fact you will not
%be able to choose a lag. It will be only a guess. The second problem is
%that these parameters were set to plot a curve(our semivariogram), 
%which will be used in the krigging process. This curve must be smooth,
%continue and must also represent the real behaviour of the semivariance.
%There are some typical models for the semivariogram adjustment. In this
%programm we will use the Gaussian Model, the   Spherical Model and the
%Exponential Model. For each model the variables are the sill and the
%range. We have already calculated the sill value and then we will obtain
%the range from these curves.
%In order to do so, a vector called DSX (DiScretized X) was implemented.
%This vector stores several values based on the distance of the input
%points. 
  
[Gaussianmodel_range, r2gaussian] = gaussianrange_evaluation(DM,LAG, sill, SVV, Weight);
[Exponentialmodel_range,r2exponential] = exponentialrange_evaluation(DM,LAG, sill, SVV, Weight );
[Sphericalmodel_range, r2spherical ]= sphericalrange_evaluation(DM,LAG, sill, SVV, Weight );

%%
%% Porosity semivariogram with three models
%The three ranges were already evaluated. This next step will plot the
%three adjustment curves and also the semivariances at each lag. Note that
%the semivariances do not need to be remade, because they were already
%computed. 


%Modeling Shperical Model to plot
clear x;
x=0:0.1:max(LAG);

for i=1:length(x)
    if x(1,i)<Sphericalmodel_range
        Sphericalmodel(i,1)=sill*(1.5*(x(1,i)/Sphericalmodel_range) - 0.5*(x(1,i)/Sphericalmodel_range)^3);
    else
        Sphericalmodel(i,1)=sill;
    end
end

%Modeling Exponential Model to plot
Exponentialmodel=sill*(1-exp((-3*x)/Exponentialmodel_range));

%Modeling Gaussian Model to plot
Gaussianmodel=sill*(1-exp((-3*power(x,2))/(Gaussianmodel_range^2)));

%Plotting all models together
figure
plot(x, Exponentialmodel);
hold on;
plot(x, Sphericalmodel);
hold on;
plot(x, Gaussianmodel);
hold on;
scatter(LAG, SVV)


%Setting Plot Legends
%Setting Graph parameters
title('Semivariograms Models')
xlabel('Lag (km)')
ylabel('Semivariance (m²/s²)')

%Setting legend parameters
legend('Exponential Model','Spherical Model', 'Gaussian Model', 'Semivariance');
legend('Location', 'southeast')
legend('Boxoff')

%%
%% Plotting the chosen model and the semivariances

figure
scatter(LAG, SVV);
hold on;
plot(x, Sphericalmodel)
title('Semivariogram - Spherical Model')
xlabel('Lag (km)')
ylabel('Semivariance (m²/s²)')

%%
%% Creating grid before kriging


griddist=100; %This number determines, indirectly, how many points will be
%interpolated, i.e., the greater this number is the finer your grid will
%be. 
DSX=linspace(0, max(A(:,1)), ((max(A(:,1))/griddist)+1));
DSY=linspace(0, max(A(:,2)),((max(A(:,2))/griddist)+1));

%This loop creates a grid by combining the discretized values of DSX and
%DSY
a=1;
for i=1:(max(A(:,1))/griddist +1)
    for j=1:(max(A(:,2))/griddist +1)
        DS(a,1)=DSX(1,i);
        DS(a,2)=DSY(1,j);
        a=a+1;
    end
end

%% Simple Kriging
[Z_sk, Z_std_sk] = kriging_simple(sill, Sphericalmodel_range, 'spherical', A, DS, DM);

%% Important Checkpoint!

%The program has just finished the 'Simple Kriging' process. It is time now
%to evaluate the ordinary kriging. But, before that, it is necessary to
%understand the differences between both. The simple kriging lies on the
%idea that the mean of the kriging variable is constant overall. The
%ordinary kriging, however, has a slightly different assumption, which is
%the idea that the variable has locals means. In other words, the mean is
%not exactlly the same over all the domain. 

%% Ordinary Kriging

[Z_ok, Z_std_ok] = kriging_ordinary(sill, Sphericalmodel_range, 'spherical', A, DS, DM);


%% 
%% Setting Plot Parameters - Variable after Simple Kriging Process

% clear x y
% x = DS(:,1);
% y = DS(:,2);

%The creation of X and Y intends to create a finer grid in order to get a 
%smoother surface at the end of the process. Still thinking of a smooth 
%surface, it is worth saying that the chosen interpolation method was the
%cubic one.
resolution=100;
range_x = min(DS(:,1)):resolution:max(DS(:,1));
range_y = min(DS(:,2)):resolution:max(DS(:,2));
[X,Y] = meshgrid(range_x,range_y);
Z = griddata(DS(:,1),DS(:,2),Z_sk,X,Y,'cubic');

figure
% surf(X,Y,Z)
imagesc(range_x,range_y,Z)
shading interp 
colorbar
colormap(jet)
set(gca,'YDir','normal')
set(gca,'XDir','normal')
caxis([min(min(Z_sk(:)),min(Z_ok(:))) max(max(Z_sk(:)),max(Z_ok(:)))]);
%The caxis command ensures that the colobar limits will be the standardize
%thecolobar limits.

%Setting legends parameters
title('Estimated Wind Speed (m/s) by Simple Kriging')
xlabel('')
ylabel('')

%% 
%% Setting Plot Parameters - Variable Standard Deviation after Simple Kriging Process

Z = griddata(DS(:,1),DS(:,2),Z_std_sk,X,Y,'cubic');

figure
% surf(X,Y,Z)
imagesc(range_x,range_y,Z)
shading interp 
colorbar
colormap(jet)
set(gca,'YDir','normal')
set(gca,'XDir','normal')
caxis([min(min(Z_std_sk(:)),min(Z_std_ok(:))) max(max(Z_std_sk(:)),max(Z_std_ok(:)))]);
%The caxis command ensures that the colobar limits will be the standardize
%thecolobar limits.

%Setting legends parameters
title('Simple Kriging Standard Deviation')
xlabel('')
ylabel('')

%% 
%% Setting Plot Parameters - Variable after Ordinary Kriging Process

%The creation of X and Y intends to create a finer grid in order to get a 
%smoother surface at the end of the process. Still thinking of a smooth 
%surface, it is worth saying that the chosen interpolation method was the
%cubic one.
Z = griddata(DS(:,1),DS(:,2),Z_ok,X,Y,'cubic');

figure
% surf(X,Y,Z)
imagesc(range_x,range_y,Z)
shading interp 
colorbar
colormap(jet)
set(gca,'YDir','normal')
set(gca,'XDir','normal')
caxis([min(min(Z_ok(:)),min(Z_ok(:))) max(max(Z_ok(:)),max(Z_ok(:)))]);
%The caxis command ensures that the colobar limits will be the standardize
%thecolobar limits.

%Setting legends parameters
title('Estimated Wind Speed (m/s) by Ordinary Kriging')
xlabel('')
ylabel('')

%% 
%% Setting Plot Parameters - Variable Standard Deviation after Ordinary Kriging Process

Z = griddata(DS(:,1),DS(:,2),Z_std_ok,X,Y,'cubic');

figure
% surf(X,Y,Z)
imagesc(range_x,range_y,Z)
shading interp 
colorbar
colormap(jet)
set(gca,'YDir','normal')
set(gca,'XDir','normal')
caxis([min(min(Z_std_sk(:)),min(Z_std_ok(:))) max(max(Z_std_sk(:)),max(Z_std_ok(:)))]);
%The caxis command ensures that the colobar limits will be the standardize
%thecolobar limits.

%Setting legends parameters
title('Ordinary Kriging Standard Deviation')
xlabel('')
ylabel('')
