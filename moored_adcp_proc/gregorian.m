function [yy,mm,dd,hh] = gregorian(jd)
%GREGORIAN Convert Julian day to Gregorian date.
%  [YY,MM,DD,HH] = GREGORIAN(JD) converts decimal Julian days to Gregorian
%  dates using the astronomical convension, but with time zero starting
%  at midnight instead of noon.  In this convention, Julian day 2440000
%  begins at 0000 hours, May 23, 1968. With one output argument a four
%  column time matrix is returned.
%
%  See also JULIAN.

%  Christian Mertens, IfM Kiel
%  adapted from the FORTRAN routine `kdate' by L. Masannek, M. Hirschberg,
%  and J. Holtorff (1980)
%  $Revision: 1.0 $ $Date: 1997/07/13 19:06:48 $

hh = rem(jd,1)*24;
jd = fix(jd-2385859);
yy = fix((4*jd - 1)/1461);
dd = fix(4*jd - 1461*yy - 1);
yy = fix(yy + 1820);
dd = fix((dd + 4)/4);
mm = fix((5*dd - 3)/153);
dd = fix(5*dd - 153*mm - 3);
dd = fix((dd + 5)/5);

mm = mm + 3;
i = (mm > 12);
yy(i) = yy(i) + 1;
mm(i) = mm(i) - 12;

if nargout <= 1
  yy = [yy mm dd hh];
end

