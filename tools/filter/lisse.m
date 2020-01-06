function [yl,yy] = lisse(y0,tab,typ,dk,fillval)
%
% Yl = lisse(Y,TAB) is a 1-D running window filter. Y = Y(i) is the signal to be filtered, 
% TAB is a vector containing the non-zero weights defining the shape of the window [i-n i+n],  
% where n*2+1 = length(TAB), n being the half-width of the window, and length(TAB) must be odd.  
% TAB allows to define manually the window size and shape. If the weights provided in TAB do not sum to 1, 
% they are normalized so that the energy of the signal is conserved. 
% 
% Yl = lisse(Y,TAB,typ) is equivalent to Yl = lisse(Y,TAB) if typ = 'mean'. 
% Using typ = 'min'('max') allows to select the minimum (maximum) value over the running window 
% rather than apply the weights.  
%
% For the n first and last points, the signal is only partially defined over the window. 
% Yl = lisse(Y,TAB,typ,dk,fillval) permits to provide a fill-value to pad the signal 
% over the window for these first and last points. 
% If fillval is the string 'extrapol', Y(1) or Y(end) are used for padding. 
% If fillval is a scalar, it is used for padding anywhere. 
% If fillval is a pair of scalars, the first value is used for padding at the beginning, 
% the second one at the end.
%
% If dk differs from 0, the signal is further shifted forward or backward and padded 
% with the corresponding fill-value, so that 
%    Yl(k) = Yl(k+dk)
% 
y = y0(:)'; 
if sum(tab)~=1, tab = tab/sum(tab); end
if nargin<5, 
  fillval_beg = NaN; 
  fillval_end = NaN; 
else
  if isstr(fillval),
    if strcmp(fillval,'extrapol'), 
      fillval_beg = y(1); fillval_end = y(end); 
    end
  else
    switch length(fillval)
     case 1, fillval_beg = fillval; fillval_end = fillval; 
     case 2, fillval_beg = fillval(1); fillval_end = fillval(2); 
    end
  end
end
if nargin<4, dk = 0; end
if nargin<3 | isempty(typ), typ = 'mean'; end
%if strcmp(typ,'max') | strcmp(typ,'min'), tab = ones(size(tab)); end
%
in = floor((length(tab)+1)/2); halfwidth = in-1; 

if size(y0,1)==1, % y0 is a vector
  %
  yy = zeros(length(tab),length(y));
  %
  switch typ
   case 'mean',
    yy(in,:) = tab(in)*y; 
    for k = 1:in-1
      decneg = in-k; 
      yy(k,:) = tab(k)*[fillval_beg*ones(1,decneg) y(1:end-decneg)]; 
    end
    for k = in+1:length(tab)
      decpos = k-in; 
      yy(k,:) = tab(k)*[y(1+decpos:end) fillval_end*ones(1,decpos)]; 
    end
    yl = nansum(yy); 
    %
   case 'mean2',
    yy(in,:) = tab(in)*y; 
    for k = 1:in-1
      decneg = in-k; 
      yy(k,:) = tab(k)*[fillval_beg*ones(1,decneg) y(1:end-decneg)].^2; 
    end
    for k = in+1:length(tab)
      decpos = k-in; 
      yy(k,:) = tab(k)*[y(1+decpos:end) fillval_end*ones(1,decpos)].^2; 
    end
    yl = sqrt(nansum(yy)); 
    %
   case 'median',
    yy(in,:) = y; 
    for k = 1:in-1
      decneg = in-k; 
      yy(k,:) = [fillval_beg*ones(1,decneg) y(1:end-decneg)]; 
    end
    for k = in+1:length(tab)
      decpos = k-in; 
      yy(k,:) = [y(1+decpos:end) fillval_end*ones(1,decpos)]; 
    end
    yl = nanmedian(yy); 
    %
   case 'max',
    yy(in,:) = y; 
    for k = 1:in-1
      decneg = in-k; 
      yy(k,:) = [fillval_beg*ones(1,decneg) y(1:end-decneg)]; 
    end
    for k = in+1:length(tab)
      decpos = k-in; 
      yy(k,:) = [y(1+decpos:end) fillval_end*ones(1,decpos)]; 
    end
    yl = nanmax(yy); 
    %
   case 'min',
    yy(in,:) = y;
    for k = 1:in-1
      decneg = in-k;
      yy(k,:) = [fillval_beg*ones(1,decneg) y(1:end-decneg)];
    end
    for k = in+1:length(tab)
      decpos = k-in;
      yy(k,:) = [y(1+decpos:end) fillval_end*ones(1,decpos)];
    end
    yl = nanmin(yy);
    %
  end
else          % y0 is a matrice
  'coucou'
end
%
y0 = yl; 
if dk == 0, yl = y0; 
elseif dk < 0, 
%keyboard
%pause
  yl = [fillval_beg*ones(1,-dk) y0(1:end+dk)]; 
elseif dk > 0, 
  yl = [y0(1+dk:end) fillval_end*ones(1,dk)]; 
end

%keyboard
%pause

return
