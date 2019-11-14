function lambda = dist_ang(lata,longa,latb,longb)
%key: distance between two points, in degrees
%synopsis :  lamdba = distance(lata,longa,latb,longb)
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
% Outputs: lamdba - distance between the two points in m
%
%uses :
%side effects :
%
%author : A.Ganachaud, Apr 95
%
%see also :
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   if any(abs(lata)>90 | abs(latb)>90)
      error('distance routine: lat > 90')
   end

   dlon=sublong(longb,longa);

   d2r=2*pi/360;
   lambda=(1/d2r)*acos( sin(d2r*lata).*sin(d2r*latb)+ ...
	cos(d2r*lata).*cos(d2r*latb).*cos(d2r*dlon) );

