function diflon = sublong( lon2, lon1 )
%key: difference between the two longitudes, according to 360deg repeat
%synopsis : diflon = sublong( lon2, lon1 )
%
%description : 
%
% diflon=lon2-lon1 ( modulo 360 )
%
%
%uses :
%
%side effects :
%
%author : A.Ganachaud (ganacho@gulf.mit.edu) , May 95
%
%see also : scan_longitude.m
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

diflon = lon2-lon1;

ii = find( diflon > 180 );
diflon(ii) = diflon(ii)-360;

jj = find( diflon < -180 );
diflon(jj) = diflon(jj)+360;
