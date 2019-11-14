%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Script:	dist.m                         Date: October 12th 1989
% Author:       A. Macdonald
% Purpose:      To determine the distance between two points
%               specified as latitude and longitude
%
% Inputs:  lata   - latitude of first point in degrees
%          longa  - longitude of the first point in degrees  
%          latb   - latitude of second point in degrees
%          longb  - longitude of the second point in degrees  
%
% Outputs: d      - distance between the two points in km
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
function d = dist(lata,longa,latb,longb)

   dlong=abs(longa-longb);

   [m,n]=size(dlong);
   for i=1:m
      if (dlong(i) > 180.) 
          dlong(i)=360.-dlong(i);
      end
   end

   dlat=lata-latb;
   
   rlata=abs(1.7453293e-2*lata);
   rlatb=abs(1.7453293e-2*latb);
   d=111.194929*sqrt(dlat.^2+(cos((rlatb+rlata)/2.).*dlong).^2);


   for i=1:m
     if((dlat(i) == 0) & (dlong(i) == 0))
        d(i)=0.0;
      end
   end
