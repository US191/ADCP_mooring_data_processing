function [par,info_rejet_med] = clean_median(param,window_lgth,nb_std,delta_lim,iter,fillval)
% par = clean_median(param,window_length,nb_std,delta_lim,iter,fillval);
% param = vector of input values
% window_lgth: number of points in the windows used to define a mean and an std_dev
% nb_std: number of std dev beyond which data are rejected
% delta_lim = [delta_min delta_max] : interval of acceptable std_dev in window_lgth 
%             (NaN possible if any value is accepted)
% iter = number of iterations of the cleaning
% fillval = value that will replace the value of the rejected data
%
% remarks: 
%  - nb_std can be a scalar or a vector with the size of 1 x iter
%  - delta_min is important in case the signal is constant over window_lgth (when std_dev=0)
%  - delta_max id also useful when the signal is very bad
%
% examples:
%    pressure = clean_median(pressure,20,[3 2.8],[1.5 10],2,-9.990e-29);
%    temp1 = clean_median(temp1,12,3,[0.05 0.4],2,-9.990e-29);
%    cond1 = clean_median(cond1,10,2.8,[0.01 0.4],3,-9.990e-29);
%    oxy1 = clean_median(oxy1,10,2.8,[0.01 0.4],3,-9.990e-29);
%
% EX: P = [1:6 40 8:20]
%  clean_median(P,5,3,[1.5 10],2,-100) : On travaille sur les donnees par paquet de 5 (=longueur de fenetre)
%  ie les donnees [1:5], puis [6 40 8 9 10], puis [11:15] puis [16:20].
% Pour chaque groupe de donnnees, on calcule la mediane et l'ecart-type sur
% ce paquet de donnees auquel on ajoute npts donnees de part et autre; npts
% etant la longueur de la demie-fenetre, soit ici npts = round(5/2) = 3.
% Ainsi, dans ce cas, la mediane et l'ecart-type sont calcules sur 11 pts
% au total (5(longueur de fenetre) + 2*npts = 5 + 2*3 = 11).
% Chacun des 5 points de la fenetre est ensuite compare avec cette mediane et cet ecart-type
% et est garde ou non.
%
% ??/??/2010 : P. Lherminier : Creation
%       2017 : P. Rousselot  : Modification valeurs absolues -> par(abs(par)<abs(med-dmax) | abs(par)>abs(med+dmax))

if length(nb_std) == 1,
    nb_std = repmat(nb_std,[iter,1]);
end
perm_param = 0;

size_param = size(param);
if size_param(1) == 1,
    param = param(:);
    perm_param = 1;
end
deja_fillval = find(param==fillval); % On repere les donnees deja a badval.
param(param==fillval) = NaN;
nparam = length(param);
nwin = floor(nparam/window_lgth)+1;
nlast = nwin*window_lgth-nparam;

par=[param;medianoutnan(param(end-nlast+1:end))*ones(nlast,1)];
par = reshape(par,[window_lgth nwin]);
% npts data are added before and after the window to calculate the std_dev (so that 
% the extremes points of the segments are not eliminated in strong gradients)
npts = round(window_lgth/2);

for iiter = 1:iter,
    parstd = [[nan(npts,1) par(end-npts+1:end,1:end-1)];par;[par(1:npts,2:end) nan(npts,1)]];
    [stdm,med] = stdmedian(parstd,0,1);
    med = med(npts+1:end-npts,:);
    dmax = repmat(min(max(nb_std(iiter)*stdm,delta_lim(1)),delta_lim(2)),[window_lgth,1]);
    par(abs(med-par)>=dmax) = NaN;            %%%%%%%%%%%%%%%%%%%%%%%
end
par = par(:); med = med(:); dmax = dmax(:);
par = par(1:nparam); med = med(1:nparam); dmax = dmax(1:nparam);

% Calcule le % de rejet apres le test d'ecart a le mediane.
nouveau_fillval = find(isnan(par));
info_rejet_med=(length(nouveau_fillval)-length(deja_fillval))*100/nparam;

par(isnan(par)) = fillval;


if perm_param == 1,
    par = par';
end

%clf
%plot(param,'b.-');
%hold on;grid on
%plot(par,'r.');
%plot([med-dmax med+dmax],'c-');

function [y,med] = stdmedian(x,flag,dim)
%STD    Standard deviation from the median and not the mean.
%   For vectors, STD(X) returns a pseudo standard deviation. For matrices,
%   STD(X) is a row vector containing the standard deviation of each
%   column.  For N-D arrays, STD(X) is the standard deviation of the
%   elements along the first non-singleton dimension of X.
%
%   STD(X) normalizes by (N-1) where N is the sequence length.  This
%   makes STD(X).^2 the best unbiased estimate of the variance if X
%   is a sample from a normal distribution.
%
%   STD(X,1) normalizes by N and produces the second moment of the
%   sample about its mean.  STD(X,0) is the same as STD(X).
%
%   STD(X,FLAG,DIM) takes the standard deviation along the dimension
%   DIM of X.  When FLAG=0 STD normalizes by (N-1), otherwise STD
%   normalizes by N.
%
%   [Y,M] = STD(X) returns the pseudo std Y and the median M
%
%   Example: If X = [4 -2 1
%                    9  5 7]
%     then std(X,0,1) is [ 3.5355 4.9497 4.2426] and std(X,0,2) is [3.0
%                                                                   2.0]
%   See also COV, MEAN, MEDIAN, CORRCOEF.

%   J.N. Little 4-21-85
%   Revised 5-9-88 JNL, 3-11-94 BAJ, 5-26-95 dlc, 5-29-96 CMT.
%   Copyright (c) 1984-98 by The MathWorks, Inc.
%   $Revision: 5.17 $  $Date: 1997/11/21 23:24:08 $

if nargin<2, flag = 0; end
if nargin<3,
  dim = find(size(x)~=1, 1 );
  if isempty(dim), dim = 1; end
end

% Avoid divide by zero.
if size(x,dim)==1, y = zeros(size(x)); return, end

tile = ones(1,max(ndims(x),dim));
tile(dim) = size(x,dim);

%xc = x - repmat(sum(x,dim)/size(x,dim),tile);  % Remove mean
med = repmat(medianoutnan(x,dim),tile);
xc = x - med;  % Remove median
xcnan=(isnan(xc));
coef=sum(~xcnan,dim);
coef(coef==0)=NaN;
xc(xcnan)=0;
if flag,
  y = sqrt(sum(conj(xc).*xc,dim)./coef);
else
  coef(coef==1)=2; %then y=0
  y = sqrt(sum(conj(xc).*xc,dim)./(coef-1));
end

function y = medianoutnan(x,dim)
%MEDIAN Median value.
%   For vectors, MEDIANOUTNAN(X) is the median value of the finite
%   elements in X. For matrices, MEDIANOUTNAN(X) is a row vector
%   containing the median value of each column.  For N-D arrays,
%   MEDIANOUTNAN(X) is the median value of the elements along the
%   first non-singleton dimension of X.
%
%   MEDIANOUTNAN(X,DIM) takes the median along the dimension DIM of X.
%
%   Example: If X = [0 1 2   NaN
%                    3 4 NaN NaN]
%
%   then medianoutnan(X,1) is [1.5 2.5 2 NaN]
%   and medianoutnan(X,2) is [1
%                             3.5]
%
%   See also MEANOUTNAN, STDOUTNAN, MIN, MAX, COV.

%   P. Lherminier, 21/03/2002. From modifications of median.m

if nargin==1,
  dim = min(find(size(x)~=1));
  if isempty(dim), dim = 1; end
end
if isempty(x), y = []; return, end

siz = [size(x) ones(1,dim-ndims(x))];
n = size(x,dim);

% Permute and reshape so that DIM becomes the row dimension of a 2-D array
perm = [dim:max(length(size(x)),dim) 1:dim-1];
x = reshape(permute(x,perm),n,prod(siz)/n);

% Sort along first dimension
x = sort(x,1);
sizx=size(x);

[nnan,ncol]=find(diff(isnan(x)));
nn=nan*ones(prod(siz)/n,1);
y=nn';
nn(ncol)=nnan;
nn(~isnan(x(end,:)))=n;

ieven=(rem(nn,2)==0);           % Even number of elements along DIM
iodd=find(~ieven & ~isnan(nn)); % Odd number of elements along DIM
ieven=find(ieven==1);
x=x(:);
if ~isempty(iodd),
    y(iodd)=x(sub2ind(sizx,(nn(iodd)+1)/2,iodd));
end;
if ~isempty(ieven),
    y(ieven)=(x(sub2ind(sizx,nn(ieven)/2,ieven))+ x(sub2ind(sizx,nn(ieven)/2+1,ieven)))/2;
end;

% Permute and reshape back
siz(dim) = 1;
y = ipermute(reshape(y,siz(perm)),perm);


