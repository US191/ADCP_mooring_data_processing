function y = median(x,dim)
% function y = median(x,dim)
%
% GEOMAR SVN $Id: nmedian.m 45 2014-06-27 15:18:59Z gkrahmann@geomar.de $
%
%NMEDIAN Median value, ignoring NaN.
%   Same as MEDIAN, but NaN's are ignored.

%   Copyright (c) 1997 by Toby Driscoll.
%   Adapted from MEDIAN.M, written by The MathWorks, Inc.
%   added backward compatibility	G.Krahmann, LODYC Paris
%   removed backward comp., added multi-dim	GK, Jun 2008

if nargin==1, 
  dim = min(find(size(x)~=1));
  if isempty(dim)
    dim = 1; 
  end
end
if isempty(x)
  y = []; 
  return
end

si = size(x);
x = shiftdim(x, ndims(x)-(dim-1));
si2 = size(x);
rx = x(:,:);
ry = repmat(nan,[1,size(rx,2)]);

for j =1:size(rx,2)
  ind = find(~isnan(rx(:,j)));
  if ~isempty(ind)
    ry(j) = median(rx(ind,j));
  end
end

si2(1) = 1;
y = reshape(ry,si2);
y = shiftdim(y,length(si2)-(ndims(x)-(dim-1)));
