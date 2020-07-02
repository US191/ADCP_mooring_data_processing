function dateNum = julianToDatenum(julian, reference)
% julianToDatenum -- Convert Julian Day to Matlab datenum.
%  julianToDatenum(dateNum) converts Julian day with days since
%   1950-01-01 00:00:00 to its equivalent dateNum (Matlab
%   datenum)  

% $Id: julianToDatenum.m 435 2010-01-14 10:31:43Z jgrelet $
%


if nargin < 1, help(mfilename), return, end

if nargin == 1
  origin = datenum(1950, 1, 1);
else
  origin = datenum(reference, 'yyyymmddHHMMSS');
end
result = origin + julian;

if nargout > 0
	dateNum = result;
else
	disp(result)
end

