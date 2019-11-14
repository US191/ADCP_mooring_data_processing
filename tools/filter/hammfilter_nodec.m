function uf = hammfilter_nodec(u,aver)

% run a box filter of length aver (points) over time series u
% replaces lowpass filter if toolbox license is missing
% no decimation (subsampling)
%
% -------- R. Zantopp  24 Nov 2006 

nn = length(u);
na = floor(aver/2);
% [w]=hamming(aver);
w  = hamm(aver);       % weights
uu = reshape(u,1,length(u));

for i = 1:nn-aver
    u1  = uu(i:i+aver-1);
    au  = u1'.*w;
    nau = sum(isnan(au));
    if nau<=floor(aver/2) 
      au_good  = find(~isnan(au));
      au2      = meanmiss(au)/meanmiss(w(au_good));
      u2(i+na) = au2;
    else
      u2(i+na) = NaN;
    end;
end;
u2(1:na)     = NaN;
u2(nn-na:nn) = NaN;
uf           = u2;
