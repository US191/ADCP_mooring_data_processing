function c = julian(yy,mm,dd,hh,m,s)

% JULIAN  Converts Gregorian calendar dates to Julian day number
%
% Although the formal definition holds that Julian days start 
%  and end at noon, here Julian days start and end at midnight.
%
% In this convention, Julian day 2440000 began at 0000 hours, May 23, 1968.
%
%  Converts Date to decimal Julian Day:
%
%    3 .. 6 Inputs:
%
%      JulianDay = JULIAN( YY , MM , DD , [hh] , [mm] , [ss] )
%
%    single Input, 3 .. 6 Columns:
%
%      JulianDay = JULIAN([ YY MM DD [hh] [mm] [ss] ]);
%
%  Converts Julian Day to Date, single Input with 1 Column:
%
%    [YY MM DD hh mm ss] = JULIAN( JulianDay )
%
%     ************************************************************
%
%        YY.... year (e.g., 1979) component
%        MM.... month (1-12) component
%        DD.... day (1-31) component of Gregorian date
%        hh.... decimal hours (assumed 0 if absent)
%        mm.... minutes
%        ss.... seconds
%
%        JulianDay  decimal Julian Day
%
%     ************************************************************
%
% see also: DATENUM, DATEVEC
%
%     recoded for MATLAB  by Rich Signell, 5-15-91
%     improved by CBegler 02/2005
%

Nin = nargin;

siz = [ 1 3 4 5 6 ];

if ~any( Nin == siz )
    error('Invalid Number of Inputs.');
end

if Nin == 1

   if isempty(yy)
      c = [];
      return
   end

   s2 = size(yy,2);

   if ~any( s2 == siz )
       error('Invalid Size of Input.');
   end

   %------------------------------------------------------
   % Day --> Date
   %------------------------------------------------------
   if s2 == 1
      c = datetime( yy + datenum(1968,05,23) - 2440000 );
      return
   end
   %------------------------------------------------------

   yy = cat( 2 , y , zeros(size(yy,1),6-size(yy,2)) );

   hh = yy(:,4) + yy(:,5)/60 + yy(:,6)/3600;
   dd = yy(:,3);
   mm = yy(:,2);
   yy = yy(:,1);

else

   if Nin < 4, hh = 0; end 
   if Nin < 5, m  = 0; end 
   if Nin < 6, s  = 0; end 

   hh = hh + m/60 + s/3600;

   sy = size(yy); py = prod(sy);
   sm = size(mm); pm = prod(sm);
   sd = size(dd); pd = prod(sd);
   sh = size(hh); ph = prod(sh);

   if ~( ( isequal(sy,sm) | ( py == 1 ) | ( pm == 1 ) ) & ...
         ( isequal(sy,sd) | ( py == 1 ) | ( pd == 1 ) ) & ...
         ( isequal(sy,sh) | ( py == 1 ) | ( ph == 1 ) ) & ...
         ( isequal(sm,sd) | ( pm == 1 ) | ( pd == 1 ) ) & ...
         ( isequal(sm,sh) | ( pm == 1 ) | ( ph == 1 ) ) & ...
         ( isequal(sd,sh) | ( pd == 1 ) | ( ph == 1 ) )       )
       error('Matrix Dimensions must be agree.');
   end

end

mo=mm+9;
yr=yy-1;

ii = ( mm > 2 );
if any(ii)
       ii  = find(ii);
    mo(ii) = mm(ii) - 3;
    yr(ii) = yy(ii);
end
 
c = floor(yr/100);

yr = yr - c*100;
      
c = floor((146097*c)/4) + floor((1461*yr)/4) + ...
    floor( ( 153*mo + 2 ) / 5 ) + dd + 1721119;

%     If you want julian days to start and end at noon, 
%     replace the following line with:
%     j=j+(h-12)/24;
 
c = c + hh/24;

%********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function d = datetime(t);

% DATETIME  returns correct DateTime
%
% Takes care on accuraccy-problems of DATEVEC, 
%    which returns seconds == 59.99999999999272
%

d      = datevec(t);

d(:,6) = round(d(:,6));  % Round seconds

dd = d(:,3);  % Original DayNumber

quot = [ 60 60 24 ]; %  [ ss-->mm  mm-->hh  hh-->dd ]
ind  = [ 6  5  4  ];

for ii = 1 : 3

    p = fix( d(:,ind(ii)) / quot(ii) );

 d(:,ind(ii)-0) = d(:,ind(ii)-0) - p * quot(ii);
 d(:,ind(ii)-1) = d(:,ind(ii)-1) + p;
  
end

% Check if DayNumber has changed

ii = find( d(:,3) > dd );

if isempty(ii)
   d(:,1:3) = d(:,1:3) + all( d(:,1:3) == 0  , 2 ) * [ -1 12 31 ];
   return
end

% New Date

[d(ii,1),d(ii,2),d(ii,3)] = datevec( datenum(d(ii,1),d(ii,2),d(ii,3)) );


d(:,1:3) = d(:,1:3) + all( d(:,1:3) == 0  , 2 ) * [ -1 12 31 ];
