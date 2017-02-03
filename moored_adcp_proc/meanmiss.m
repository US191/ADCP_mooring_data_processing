function y = meanmiss(x,dim)
%MEANMISS Average or mean value with missing data.
%  Y = MEANMISS(X)
%  Y = MEANMISS(X,DIM)


if nargin ==1
  dim = min(find(size(x) ~= 1));
  if isempty(dim)
    dim = 1;
  end
end

y = zeros(size(x));
valid = ~isnan(x);
y(valid) = x(valid);
n = sum(valid,dim);
n(n==0) = NaN;
y = sum(y,dim)./n;

