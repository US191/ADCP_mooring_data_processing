function d = distance(lata,longa,latb,longb)
%key: distance between two points, in METERS
%synopsis : d = distance(lata,longa,latb,longb)
%
%description : 
%
% Purpose:      To determine the distance between two points
%               specified as latitude and longitude
%
% Inputs:  lata   - latitude of first point in degrees
%          longa  - longitude of the first point in degrees  
%          latb   - latitude of second point in degrees
%          longb  - longitude of the second point in degrees  
%
% Outputs: d      - distance between the two points in m
%
%uses : the latitudes must not be too different
%
%side effects :
%
%author : A. Macdonald -  A.Ganachaud, Apr 95
%
%see also :
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   if any(abs(lata)>90 | abs(latb)>90)
      error('distance routine: lat > 90')
   end
   dlong=abs(longa-longb);

   dbad=find(dlong >= 180.); 
   dlong(dbad)=360-dlong(dbad);

   dlat=lata-latb;
   
   rlata=abs(1.7453293e-2*lata);
   rlatb=abs(1.7453293e-2*latb);
   d=111194.929*sqrt(dlat.^2+(cos((rlatb+rlata)/2.).*dlong).^2);


   dzer=find((dlat == 0) & (dlong == 0));
   d(dzer)=zeros(size(dzer));