function [DM] = distance_measurement(A, type)

DM=zeros(length(A'));

if strcmpi(type, 'geodesic') ==1 
    
%Once we are dealing with global coordinates in a very large space it 
%will be needed to evaluate the geodesic distance. Matlab has a very
%helpfull feature which allows to set a reference ellipsoid and evaluate the
%distance once the latitude and longitude of each point were given
    
    Earth = referenceEllipsoid('Geodetic Reference System 1980', 'km');
    for i = 1:length(A')
       for j = 1:length(A');
          DM(i,j) = distance(A(i,1),A(i,2),A(j,1),A(j,2), Earth); 
       end
    end

elseif strcmpi(type, 'cartesian') == 1
    
    AUX = A;
    AUX(:,3) = [];
    DM = squareform(pdist(AUX));
    
else
    error;
        
DM = DM;
    
end