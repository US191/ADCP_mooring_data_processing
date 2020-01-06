function varargout = t_signal(fcn,varargin);

% T_SIGNAL  PSD and CSD from Signal Processing toolbox
%
% T_SIGNAL( Fcn , varargin )
%

%---------------------------------------------------------------------
% SubFunctions from Signal and Matlab-Toolboxes
%---------------------------------------------------------------------
%
% PSD        Power Spectral Density estimate.
% CSD        Cross Spectral Density estimate.
%
% PSDCHK     Helper function for PSD, CSD, COHERE, and TFE.
%
% CHI2CONF   Confidence interval using inverse of chi-square cdf.
% -CHI2INV   Inverse of the chi-square cumulative distribution function (cdf).
% -GAMINV    Inverse of the gamma cumulative distribution function (cdf).
% -NORMINV   Inverse of the normal cumulative distribution function (cdf).
% -GAMCDF    Gamma cumulative distribution function.
% -GAMPDF    Gamma probability density function.
% -DISTCHCK  Checks the argument list for the probability functions.
%
% HAMMING       Hamming window.
%
% GENCOSWIN     Returns one of the generalized cosine windows.
% -SYM_WINDOW   Symmetric generalized cosine window.
% -CALC_WINDOW  Calculate the generalized cosine window samples.
%
% HANNING       Hanning window.
% -SYM_HANNING  Symmetric Hanning window. 
% -CALC_HANNING Calculates Hanning window samples.
%
% CHECK_ORDER   Checks the order passed to the window functions.
%

varargout = cell(1,nargout);

[varargout{:}] = feval(['t_' fcn],varargin{:});

%*****************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [Pxx, Pxxc, f] = t_psd(varargin)
%PSD Power Spectral Density estimate.
%   PSD has been replaced by SPECTRUM.WELCH.  PSD still works but may be
%   removed in the future. Use SPECTRUM.WELCH (or its functional form
%   PWELCH) instead. Type help SPECTRUM/WELCH for details.
%
%   See also SPECTRUM/PSD, SPECTRUM/MSSPECTRUM, SPECTRUM/PERIODOGRAM.

%   Author(s): T. Krauss, 3-26-93
%   Copyright 1988-2004 The MathWorks, Inc.
%   $Revision: 1.12.4.3 $  $Date: 2004/10/18 21:09:08 $

%   NOTE 1: To express the result of PSD, Pxx, in units of
%           Power per Hertz multiply Pxx by 1/Fs [1].
%
%   NOTE 2: The Power Spectral Density of a continuous-time signal,
%           Pss (watts/Hz), is proportional to the Power Spectral 
%           Density of the sampled discrete-time signal, Pxx, by Ts
%           (sampling period). [2] 
%       
%               Pss(w/Ts) = Pxx(w)*Ts,    |w| < pi; where w = 2*pi*f*Ts

%   References:
%     [1] Petre Stoica and Randolph Moses, Introduction To Spectral
%         Analysis, Prentice hall, 1997, pg, 15
%     [2] A.V. Oppenheim and R.W. Schafer, Discrete-Time Signal
%         Processing, Prentice-Hall, 1989, pg. 731
%     [3] A.V. Oppenheim and R.W. Schafer, Digital Signal
%         Processing, Prentice-Hall, 1975, pg. 556

error(nargchk(1,7,nargin))
x = varargin{1};
[msg,nfft,Fs,window,noverlap,p,dflag]=psdchk(varargin(2:end),x);
error(msg)

% compute PSD
window = window(:);
n = length(x);		    % Number of data points
nwind = length(window); % length of window
if n < nwind            % zero-pad x if it has length less than the window length
    x(nwind)=0;  n=nwind;
end
% Make sure x is a column vector; do this AFTER the zero-padding 
% in case x is a scalar.
x = x(:);		

k = fix((n-noverlap)/(nwind-noverlap));	% Number of windows
                    					% (k = fix(n/nwind) for noverlap=0)
%   if 0
%       disp(sprintf('   x        = (length %g)',length(x)))
%       disp(sprintf('   y        = (length %g)',length(y)))
%       disp(sprintf('   nfft     = %g',nfft))
%       disp(sprintf('   Fs       = %g',Fs))
%       disp(sprintf('   window   = (length %g)',length(window)))
%       disp(sprintf('   noverlap = %g',noverlap))
%       if ~isempty(p)
%           disp(sprintf('   p        = %g',p))
%       else
%           disp('   p        = undefined')
%       end
%       disp(sprintf('   dflag    = ''%s''',dflag))
%       disp('   --------')
%       disp(sprintf('   k        = %g',k))
%   end

index = 1:nwind;
KMU = k*norm(window)^2;	% Normalizing scale factor ==> asymptotically unbiased
% KMU = k*sum(window)^2;% alt. Nrmlzng scale factor ==> peaks are about right

Spec = zeros(nfft,1); % Spec2 = zeros(nfft,1);
for i=1:k
    if strcmp(dflag,'none')
        xw = window.*(x(index));
    elseif strcmp(dflag,'linear')
        xw = window.*detrend(x(index));
    else
        xw = window.*detrend(x(index),'constant');
    end
    index = index + (nwind - noverlap);
    Xx = abs(fft(xw,nfft)).^2;
    Spec = Spec + Xx;
%     Spec2 = Spec2 + abs(Xx).^2;
end

% Select first half
if ~any(any(imag(x)~=0)),   % if x is not complex
    if rem(nfft,2),    % nfft odd
        select = (1:(nfft+1)/2)';
    else
        select = (1:nfft/2+1)';
    end
    Spec = Spec(select);
%     Spec2 = Spec2(select);
%    Spec = 4*Spec(select);     % double the signal content - essentially
% folding over the negative frequencies onto the positive and adding.
%    Spec2 = 16*Spec2(select);
else
    select = (1:nfft)';
end
freq_vector = (select - 1)*Fs/nfft;

% find confidence interval if needed
if (nargout == 3) || ((nargout == 0) && ~isempty(p)),
    if isempty(p),
        p = .95;    % default
    end
    % Confidence interval from Kay, p. 76, eqn 4.16:
    % (first column is lower edge of conf int., 2nd col is upper edge)
    confid = Spec*chi2conf(p,k)/KMU;

    if noverlap > 0
        disp('Warning: confidence intervals inaccurate for NOVERLAP > 0.')
    end
end

Spec = Spec*(1/KMU);   % normalize

% set up output parameters
if (nargout == 3),
   Pxx = Spec;
   Pxxc = confid;
   f = freq_vector;
elseif (nargout == 2),
   Pxx = Spec;
   Pxxc = freq_vector;
elseif (nargout == 1),
   Pxx = Spec;
elseif (nargout == 0),
   if ~isempty(p),
       P = [Spec confid];
   else
       P = Spec;
   end
   newplot;
   plot(freq_vector,10*log10(abs(P))), grid on
   xlabel('Frequency'), ylabel('Power Spectrum Magnitude (dB)');
end

%*****************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [Pxy, Pxyc, f] = t_csd(varargin)
%CSD Cross Spectral Density estimate.
%   CSD has been replaced by CPSD.  CSD still works but may be removed in
%   the future. Use CPSD instead.
%
%   CPSD does not support detrending.  Please use the DETREND function if
%   you need to detrend your signal.  Type "help detrend" for details.
%
%   See also PWELCH, TFESTIMATE, MSCOHERE,
%   ETFE, SPA, and ARX in the Identification Toolbox.

%   Author(s): T. Krauss, 3-30-93
%   Copyright 1988-2004 The MathWorks, Inc.
%   $Revision: 1.8.4.4 $  $Date: 2004/12/26 22:15:35 $

error(nargchk(2,8,nargin))
x = varargin{1};
y = varargin{2};
[msg,nfft,Fs,window,noverlap,p,dflag]=psdchk(varargin(3:end),x,y);
error(msg)
    
% compute CSD
window = window(:);
n = length(x);		% Number of data points
nwind = length(window); % length of window
if n < nwind    % zero-pad x , y if length is less than the window length
    x(nwind)=0;
    y(nwind)=0;  
    n=nwind;
end
x = x(:);		% Make sure x is a column vector
y = y(:);		% Make sure y is a column vector
k = fix((n-noverlap)/(nwind-noverlap));	% Number of windows
					% (k = fix(n/nwind) for noverlap=0)
index = 1:nwind;
KMU = k*norm(window)^2;	% Normalizing scale factor ==> asymptotically unbiased
% KMU = k*sum(window)^2;% alt. Nrmlzng scale factor ==> peaks are about right

Spec = zeros(nfft,1);
for i=1:k
    if strcmp(dflag,'none')
        xw = window.*x(index);
        yw = window.*y(index);
    elseif strcmp(dflag,'linear')
        xw = window.*detrend(x(index));
        yw = window.*detrend(y(index));
    else
        xw = window.*detrend(x(index),0);
        yw = window.*detrend(y(index),0);
    end
    index = index + (nwind - noverlap);
    Xx = fft(xw,nfft);
    Yy = fft(yw,nfft);
    Xy2 = Yy.*conj(Xx);
    Spec = Spec + Xy2;
end

% Select first half
if ~any(any(imag([x y])~=0)),   % if x and y are not complex
    if rem(nfft,2),    % nfft odd
        select = 1:(nfft+1)/2;
    else
        select = 1:nfft/2+1;   % include DC AND Nyquist
    end
    Spec = Spec(select);
else
    select = 1:nfft;
end
freq_vector = (select - 1)'*Fs/nfft;

% find confidence interval if needed
if (nargout == 3) || ((nargout == 0) && ~isempty(p)),
    if isempty(p),
        p = .95;    % default
    end
    confid = Spec*chi2conf(p,k)/KMU;
end

Spec = Spec*(1/KMU);

% set up output parameters
if (nargout == 3),
   Pxy = Spec;
   Pxyc = confid;
   f = freq_vector;
elseif (nargout == 2),
   Pxy = Spec;
   Pxyc = freq_vector;
elseif (nargout == 1),
   Pxy = Spec;
elseif (nargout == 0),
   if ~isempty(p),
       P = [Spec confid];
   else
       P = Spec;
   end
   newplot;
   plot(freq_vector,10*log10(abs(P))), grid on
   xlabel('Frequency'), ylabel('Cross Spectrum Magnitude (dB)');
end

%*****************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [msg,nfft,Fs,window,noverlap,p,dflag] = psdchk(P,x,y)
%PSDCHK Helper function for PSD, CSD, COHERE, and TFE.
%   [msg,nfft,Fs,window,noverlap,p,dflag]=PSDCHK(P,x,y) takes the cell 
%   array P and uses each element as an input argument.  Assumes P has 
%   between 0 and 7 elements which are the arguments to psd, csd, cohere
%   or tfe after the x (psd) or x and y (csd, cohere, tfe) arguments.
%   y is optional; if given, it is checked to match the size of x.
%   x must be a numeric vector.
%   Outputs:
%     msg - error message, [] if no error
%     nfft - fft length
%     Fs - sampling frequency
%     window - window vector
%     noverlap - overlap of sections, in samples
%     p - confidence interval, [] if none desired
%     dflag - detrending flag, 'linear' 'mean' or 'none'

%   Author(s): T. Krauss, 10-28-93
%   Copyright 1988-2002 The MathWorks, Inc.
%       $Revision: 1.6 $  $Date: 2002/04/15 01:08:06 $

msg = [];

if length(P) == 0 
% psd(x)
    nfft = min(length(x),256);
    window = hanning(nfft);
    noverlap = 0;
    Fs = 2;
    p = [];
    dflag = 'none';
elseif length(P) == 1
% psd(x,nfft)
% psd(x,dflag)
    if isempty(P{1}),   dflag = 'none'; nfft = min(length(x),256); 
    elseif isstr(P{1}), dflag = P{1};       nfft = min(length(x),256); 
    else              dflag = 'none'; nfft = P{1};   end
    Fs = 2;
    window = hanning(nfft);
    noverlap = 0;
    p = [];
elseif length(P) == 2
% psd(x,nfft,Fs)
% psd(x,nfft,dflag)
    if isempty(P{1}), nfft = min(length(x),256); else nfft=P{1};     end
    if isempty(P{2}),   dflag = 'none'; Fs = 2;
    elseif isstr(P{2}), dflag = P{2};       Fs = 2;
    else              dflag = 'none'; Fs = P{2}; end
    window = hanning(nfft);
    noverlap = 0;
    p = [];
elseif length(P) == 3
% psd(x,nfft,Fs,window)
% psd(x,nfft,Fs,dflag)
    if isempty(P{1}), nfft = min(length(x),256); else nfft=P{1};     end
    if isempty(P{2}), Fs = 2;     else    Fs = P{2}; end
    if isstr(P{3})
        dflag = P{3};
        window = hanning(nfft);
    else
        dflag = 'none';
        window = P{3};
        if length(window) == 1, window = hanning(window); end
        if isempty(window), window = hanning(nfft); end
    end
    noverlap = 0;
    p = [];
elseif length(P) == 4
% psd(x,nfft,Fs,window,noverlap)
% psd(x,nfft,Fs,window,dflag)
    if isempty(P{1}), nfft = min(length(x),256); else nfft=P{1};     end
    if isempty(P{2}), Fs = 2;     else    Fs = P{2}; end
    window = P{3};
    if length(window) == 1, window = hanning(window); end
    if isempty(window), window = hanning(nfft); end
    if isstr(P{4})
        dflag = P{4};
        noverlap = 0;
    else
        dflag = 'none';
        if isempty(P{4}), noverlap = 0; else noverlap = P{4}; end
    end
    p = [];
elseif length(P) == 5
% psd(x,nfft,Fs,window,noverlap,p)
% psd(x,nfft,Fs,window,noverlap,dflag)
    if isempty(P{1}), nfft = min(length(x),256); else nfft=P{1};     end
    if isempty(P{2}), Fs = 2;     else    Fs = P{2}; end
    window = P{3};
    if length(window) == 1, window = hanning(window); end
    if isempty(window), window = hanning(nfft); end
    if isempty(P{4}), noverlap = 0; else noverlap = P{4}; end
    if isstr(P{5})
        dflag = P{5};
        p = [];
    else
        dflag = 'none';
        if isempty(P{5}), p = .95;    else    p = P{5}; end
    end
elseif length(P) == 6
% psd(x,nfft,Fs,window,noverlap,p,dflag)
    if isempty(P{1}), nfft = min(length(x),256); else nfft=P{1};     end
    if isempty(P{2}), Fs = 2;     else    Fs = P{2}; end
    window = P{3};
    if length(window) == 1, window = hanning(window); end
    if isempty(window), window = hanning(nfft); end
    if isempty(P{4}), noverlap = 0; else noverlap = P{4}; end
    if isempty(P{5}), p = .95;    else    p = P{5}; end
    if isstr(P{6})
        dflag = P{6};
    else
        msg = 'DFLAG parameter must be a string.'; return
    end
end

% NOW do error checking
if (nfft<length(window)), 
    msg = 'Requires window''s length to be no greater than the FFT length.';
end
if (noverlap >= length(window)),
    msg = 'Requires NOVERLAP to be strictly less than the window length.';
end
if (nfft ~= abs(round(nfft)))|(noverlap ~= abs(round(noverlap))),
    msg = 'Requires positive integer values for NFFT and NOVERLAP.';
end
if ~isempty(p),
    if (prod(size(p))>1)|(p(1,1)>1)|(p(1,1)<0),
        msg = 'Requires confidence parameter to be a scalar between 0 and 1.';
    end
end
if min(size(x))~=1 | ~isnumeric(x) | length(size(x))>2
    msg = 'Requires vector (either row or column) input.';
end
if (nargin>2) & ( (min(size(y))~=1) | ~isnumeric(y) | length(size(y))>2 )
    msg = 'Requires vector (either row or column) input.';
end
if (nargin>2) & (length(x)~=length(y))
    msg = 'Requires X and Y be the same length.';
end

dflag = lower(dflag);
if strncmp(dflag,'none',1)
      dflag = 'none';
elseif strncmp(dflag,'linear',1)
      dflag = 'linear';
elseif strncmp(dflag,'mean',1)
      dflag = 'mean';
else
    msg = 'DFLAG must be ''linear'', ''mean'', or ''none''.';
end

%*****************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function c = chi2conf(conf,k);
%CHI2CONF Confidence interval using inverse of chi-square cdf.
%   C = CHI2CONF(P,K) is the confidence interval of an unbiased power spectrum 
%   estimate made up of K independent measurements.  C is a two element
%   vector.  We are P*100% confident that the true PSD lies in the interval
%   [C(1)*X C(2)*X], where X is the PSD estimate.
%
%   Reference:
%     Stephen Kay, "Modern Spectral Analysis, Theory & Application," 
%     p. 76, eqn 4.16.

%   Copyright 1988-2002 The MathWorks, Inc.
%   $Revision: 1.6 $  $Date: 2002/04/15 01:07:39 $

if nargin < 2, 
    error('Requires two input arguments.');
end

v=2*k;
alfa = 1 - conf;
c=chi2inv([1-alfa/2 alfa/2],v);
c=v./c;

function x = chi2inv(p,v);
%CHI2INV Inverse of the chi-square cumulative distribution function (cdf).
%   X = CHI2INV(P,V)  returns the inverse of the chi-square cdf with V  
%   degrees of freedom at the values in P. The chi-square cdf with V 
%   degrees of freedom, is the gamma cdf with parameters V/2 and 2.   
%
%   The size of X is the common size of P and V. A scalar input
%   functions as a constant matrix of the same size as the other input.   

%   References:
%      [1]  M. Abramowitz and I. A. Stegun, "Handbook of Mathematical
%      Functions", Government Printing Office, 1964, 26.4.
%      [2] E. Kreyszig, "Introductory Mathematical Statistics",
%      John Wiley, 1970, section 10.2 (page 144)

%  Was: Revision: 1.2, Date: 1996/07/25 16:23:36

if nargin < 2, 
    error('Requires two input arguments.');
end

[errorcode p v] = distchck(2,p,v);

if errorcode > 0
    error('Requires non-scalar arguments to match in size.');
end

% Call the gamma inverse function. 
x = gaminv(p,v/2,2);

% Return NaN if the degrees of freedom is not a positive integer.
k = find(v < 0  |  round(v) ~= v);
if any(k)
    tmp  = NaN;
    x(k) = tmp(ones(size(k)));
end


function x = gaminv(p,a,b);
%GAMINV Inverse of the gamma cumulative distribution function (cdf).
%   X = GAMINV(P,A,B)  returns the inverse of the gamma cdf with  
%   parameters A and B, at the probabilities in P.
%
%   The size of X is the common size of the input arguments. A scalar input  
%   functions as a constant matrix of the same size as the other inputs.    
%
%   GAMINV uses Newton's method to converge to the solution.

%   References:
%      [1]  M. Abramowitz and I. A. Stegun, "Handbook of Mathematical
%      Functions", Government Printing Office, 1964, 6.5.

%   B.A. Jones 1-12-93
%   Was: Revision: 1.2, Date: 1996/07/25 16:23:36

if nargin<3, 
    b=1;
end

[errorcode p a b] = distchck(3,p,a,b);

if errorcode > 0
    error('The arguments must be the same size or be scalars.');
end

%   Initialize X to zero.
x = zeros(size(p));

k = find(p<0 | p>1 | a <= 0 | b <= 0);
if any(k),
    tmp = NaN;
    x(k) = tmp(ones(size(k)));
end

% The inverse cdf of 0 is 0, and the inverse cdf of 1 is 1.  
k0 = find(p == 0 & a > 0 & b > 0);
if any(k0),
    x(k0) = zeros(size(k0)); 
end

k1 = find(p == 1 & a > 0 & b > 0);
if any(k1), 
    tmp = Inf;
    x(k1) = tmp(ones(size(k1))); 
end

% Newton's Method
% Permit no more than count_limit interations.
count_limit = 100;
count = 0;

k = find(p > 0  &  p < 1 & a > 0 & b > 0);
pk = p(k);

% Supply a starting guess for the iteration.
%   Use a method of moments fit to the lognormal distribution. 
mn = a(k) .* b(k);
v = mn .* b(k);
temp = log(v + mn .^ 2); 
mu = 2 * log(mn) - 0.5 * temp;
sigma = -2 * log(mn) + temp;
xk = exp(norminv(pk,mu,sigma));

h = ones(size(pk)); 

% Break out of the iteration loop for three reasons:
%  1) the last update is very small (compared to x)
%  2) the last update is very small (compared to sqrt(eps))
%  3) There are more than 100 iterations. This should NEVER happen. 

while(any(abs(h) > sqrt(eps)*abs(xk))  &  max(abs(h)) > sqrt(eps)    ...
                                 & count < count_limit), 
                                 
    count = count + 1;
    h = (gamcdf(xk,a(k),b(k)) - pk) ./ gampdf(xk,a(k),b(k));
    xnew = xk - h;
    % Make sure that the current guess stays greater than zero.
    % When Newton's Method suggests steps that lead to negative guesses
    % take a step 9/10ths of the way to zero:
    ksmall = find(xnew < 0);
    if any(ksmall),
        xnew(ksmall) = xk(ksmall) / 10;
        h = xk-xnew;
    end
    xk = xnew;
end


% Store the converged value in the correct place
x(k) = xk;

if count == count_limit, 
    fprintf('\nWarning: GAMINV did not converge.\n');
    str = 'The last step was:  ';
    outstr = sprintf([str,'%13.8f'],h);
    fprintf(outstr);
end

function x = norminv(p,mu,sigma);
%NORMINV Inverse of the normal cumulative distribution function (cdf).
%   X = NORMINV(P,MU,SIGMA) finds the inverse of the normal cdf with
%   mean, MU, and standard deviation, SIGMA.
%
%   The size of X is the common size of the input arguments. A scalar input  
%   functions as a constant matrix of the same size as the other inputs.    
%
%   Default values for MU and SIGMA are 0 and 1 respectively.

%   References:
%      [1]  M. Abramowitz and I. A. Stegun, "Handbook of Mathematical
%      Functions", Government Printing Office, 1964, 7.1.1 and 26.2.2

%   Was: Revision: 1.2, Date: 1996/07/25 16:23:36

if nargin < 3, 
    sigma = 1;
end

if nargin < 2;
    mu = 0;
end

[errorcode p mu sigma] = distchck(3,p,mu,sigma);

if errorcode > 0
    error('Requires non-scalar arguments to match in size.');
end

% Allocate space for x.
x = zeros(size(p));

% Return NaN if the arguments are outside their respective limits.
k = find(sigma <= 0 | p < 0 | p > 1);
if any(k)
    tmp  = NaN;
    x(k) = tmp(ones(size(k))); 
end

% Put in the correct values when P is either 0 or 1.
k = find(p == 0);
if any(k)
    tmp  = Inf;
    x(k) = -tmp(ones(size(k)));
end

k = find(p == 1);
if any(k)
    tmp  = Inf;
    x(k) = tmp(ones(size(k))); 
end

% Compute the inverse function for the intermediate values.
k = find(p > 0  &  p < 1 & sigma > 0);
if any(k),
    x(k) = sqrt(2) * sigma(k) .* erfinv(2 * p(k) - 1) + mu(k);
end


function p = gamcdf(x,a,b);
%GAMCDF Gamma cumulative distribution function.
%   P = GAMCDF(X,A,B) returns the gamma cumulative distribution
%   function with parameters A and B at the values in X.
%
%   The size of P is the common size of the input arguments. A scalar input  
%   functions as a constant matrix of the same size as the other inputs.    
%
%   Some references refer to the gamma distribution with a single
%   parameter. This corresponds to the default of B = 1. 
%
%   GAMMAINC does computational work.

%   References:
%      [1]  L. Devroye, "Non-Uniform Random Variate Generation", 
%      Springer-Verlag, 1986. p. 401.
%      [2]  M. Abramowitz and I. A. Stegun, "Handbook of Mathematical
%      Functions", Government Printing Office, 1964, 26.1.32.

%   Was: Revision: 1.2, Date: 1996/07/25 16:23:36

if nargin < 3, 
    b = 1; 
end

if nargin < 2, 
    error('Requires at least two input arguments.'); 
end

[errorcode x a b] = distchck(3,x,a,b);

if errorcode > 0
    error('Requires non-scalar arguments to match in size.');
end

%   Return NaN if the arguments are outside their respective limits.
k1 = find(a <= 0 | b <= 0);     
if any(k1)
    tmp   = NaN;
    p(k1) = tmp(ones(size(k1)));
end

% Initialize P to zero.
p = zeros(size(x));

k = find(x > 0 & ~(a <= 0 | b <= 0));
if any(k), 
    p(k) = gammainc(x(k) ./ b(k),a(k));
end

% Make sure that round-off errors never make P greater than 1.
k = find(p > 1);
if any(k)
    p(k) = ones(size(k));
end


function y = gampdf(x,a,b)
%GAMPDF Gamma probability density function.
%   Y = GAMPDF(X,A,B) returns the gamma probability density function 
%   with parameters A and B, at the values in X.
%
%   The size of Y is the common size of the input arguments. A scalar input  
%   functions as a constant matrix of the same size as the other inputs.    
%
%   Some references refer to the gamma distribution with a single
%   parameter. This corresponds to the default of B = 1.

%   References:
%      [1]  L. Devroye, "Non-Uniform Random Variate Generation", 
%      Springer-Verlag, 1986, pages 401-402.

%   Was: Revision: 1.2, Date: 1996/07/25 16:23:36

if nargin < 3, 
    b = 1; 
end

if nargin < 2, 
    error('Requires at least two input arguments'); 
end

[errorcode x a b] = distchck(3,x,a,b);

if errorcode > 0
    error('Requires non-scalar arguments to match in size.');
end

% Initialize Y to zero.
y = zeros(size(x));

%   Return NaN if the arguments are outside their respective limits.
k1 = find(a <= 0 | b <= 0);     
if any(k1)
    tmp = NaN;
    y(k1) = tmp(ones(size(k1)));
end

k=find(x > 0 & ~(a <= 0 | b <= 0));
if any(k)
    y(k) = (a(k) - 1) .* log(x(k)) - (x(k) ./ b(k)) - gammaln(a(k)) - a(k) .* log(b(k));
    y(k) = exp(y(k));
end
k1 = find(x == 0 & a < 1);
if any(k1)
  tmp = Inf;
  y(k1) = tmp(ones(size(k1)));
end
k2 = find(x == 0 & a == 1);
if any(k2)
  y(k2) = (1./b(k2));
end

function [errorcode,out1,out2,out3,out4] = distchck(nparms,arg1,arg2,arg3,arg4)
%DISTCHCK Checks the argument list for the probability functions.

%   B.A. Jones  1-22-93
%   Was: Revision: 1.2, Date: 1996/07/25 16:23:36

errorcode = 0;

if nparms == 1
    out1 = arg1;
    return;
end
    
if nparms == 2
    [r1 c1] = size(arg1);
    [r2 c2] = size(arg2);
    scalararg1 = (prod(size(arg1)) == 1);
    scalararg2 = (prod(size(arg2)) == 1);
    if ~scalararg1 & ~scalararg2
        if r1 ~= r2 | c1 ~= c2
            errorcode = 1;
            return;         
        end     
    end
    if scalararg1
        out1 = arg1(ones(r2,1),ones(c2,1));
    else
        out1 = arg1;
    end
    if scalararg2
        out2 = arg2(ones(r1,1),ones(c1,1));
    else
        out2 = arg2;
    end
end
    
if nparms == 3
    [r1 c1] = size(arg1);
    [r2 c2] = size(arg2);
    [r3 c3] = size(arg3);
    scalararg1 = (prod(size(arg1)) == 1);
    scalararg2 = (prod(size(arg2)) == 1);
    scalararg3 = (prod(size(arg3)) == 1);

    if ~scalararg1 & ~scalararg2
        if r1 ~= r2 | c1 ~= c2
            errorcode = 1;
            return;         
        end
    end

    if ~scalararg1 & ~scalararg3
        if r1 ~= r3 | c1 ~= c3
            errorcode = 1;
            return;                 
        end
    end

    if ~scalararg3 & ~scalararg2
        if r3 ~= r2 | c3 ~= c2
            errorcode = 1;
            return;         
        end
    end

    if ~scalararg1
      out1 = arg1;
    end
    if ~scalararg2
      out2 = arg2;
    end
    if ~scalararg3
      out3 = arg3;
    end
    rows = max([r1 r2 r3]);
   columns = max([c1 c2 c3]);  
       
    if scalararg1
        out1 = arg1(ones(rows,1),ones(columns,1));
    end
   if scalararg2
        out2 = arg2(ones(rows,1),ones(columns,1));
   end
   if scalararg3
       out3 = arg3(ones(rows,1),ones(columns,1));
   end
     out4 =[];
    
end

if nparms == 4
    [r1 c1] = size(arg1);
    [r2 c2] = size(arg2);
    [r3 c3] = size(arg3);
    [r4 c4] = size(arg4);
    scalararg1 = (prod(size(arg1)) == 1);
    scalararg2 = (prod(size(arg2)) == 1);
    scalararg3 = (prod(size(arg3)) == 1);
    scalararg4 = (prod(size(arg4)) == 1);

    if ~scalararg1 & ~scalararg2
        if r1 ~= r2 | c1 ~= c2
            errorcode = 1;
            return;         
        end
    end

    if ~scalararg1 & ~scalararg3
        if r1 ~= r3 | c1 ~= c3
            errorcode = 1;
            return;                 
        end
    end

    if ~scalararg1 & ~scalararg4
        if r1 ~= r4 | c1 ~= c4
            errorcode = 1;
            return;                 
        end
    end

    if ~scalararg3 & ~scalararg2
        if r3 ~= r2 | c3 ~= c2
            errorcode = 1;
            return;         
        end
    end

    if ~scalararg4 & ~scalararg2
        if r4 ~= r2 | c4 ~= c2
            errorcode = 1;
            return;         
        end
    end

    if ~scalararg3 & ~scalararg4
        if r3 ~= r4 | c3 ~= c4
            errorcode = 1;
            return;         
        end
    end


    if ~scalararg1
       out1 = arg1;
    end
    if ~scalararg2
       out2 = arg2;
    end
    if ~scalararg3
      out3 = arg3;
    end
    if ~scalararg4
      out4 = arg4;
    end
 
   rows = max([r1 r2 r3 r4]);
   columns = max([c1 c2 c3 c4]);     
    if scalararg1
       out1 = arg1(ones(rows,1),ones(columns,1));
   end
   if scalararg2
       out2 = arg2(ones(rows,1),ones(columns,1));
   end
   if scalararg3
       out3 = arg3(ones(rows,1),ones(columns,1));
   end
   if scalararg4
       out4 = arg4(ones(rows,1),ones(columns,1));
   end
end


%*****************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function w = hamming(varargin)
%HAMMING   Hamming window.
%   HAMMING(N) returns the N-point symmetric Hamming window in a column vector.
% 
%   HAMMING(N,SFLAG) generates the N-point Hamming window using SFLAG window
%   sampling. SFLAG may be either 'symmetric' or 'periodic'. By default, a 
%   symmetric window is returned. 
%
%   See also BLACKMAN, HANN, WINDOW.

%   Copyright 1988-2002 The MathWorks, Inc.
%   $Revision: 1.14 $  $Date: 2002/11/21 15:46:43 $

% Check number of inputs
error(nargchk(1,2,nargin));

[w,msg] = gencoswin('hamming',varargin{:});
error(msg);


% [EOF] hamming.m

%*****************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [w,msg] = gencoswin(varargin)
%GENCOSWIN   Returns one of the generalized cosine windows.
%   GENCOSWIN returns the generalized cosine window specified by the 
%   first string argument. Its inputs can be
%     Window name    - a string, any of 'hamming', 'hann', 'blackman'.
%     N              - length of the window desired.
%     Sampling flag  - optional string, one of 'symmetric', 'periodic'. 

%   Copyright 1988-2002 The MathWorks, Inc.
%   $Revision: 1.7 $  $Date: 2002/04/15 01:09:14 $ 

% Parse the inputs
window = varargin{1};
n = varargin{2};
msg = '';

% Check for trivial orders
[n,w,trivialwin] = check_order(n);
if trivialwin, return, end;

% Select the sampling option
if nargin == 2, % no sampling flag specified, use default. 
    sflag = 'symmetric';
else
    sflag = lower(varargin{3});
end

% Allow partial strings for sampling options
allsflags = {'symmetric','periodic'};
sflagindex = strmatch(sflag, allsflags);
if length(sflagindex)~=1         % catch 0 or 2 matches
    msg = 'Sampling flag must be either ''symmetric'' or ''periodic''.';
    return;
else	
    sflag = allsflags{sflagindex};
end

% Evaluate the window
switch sflag
case 'periodic'
    w = sym_window(n+1,window);
    w(end) = [];
case 'symmetric'
    w = sym_window(n,window);
end

%---------------------------------------------------------------------
function w = sym_window(n,window)
%SYM_WINDOW   Symmetric generalized cosine window.
%   SYM_WINDOW Returns an exactly symmetric N point generalized cosine 
%   window by evaluating the first half and then flipping the same samples
%   over the other half.

if ~rem(n,2)
    % Even length window
    half = n/2;
    w = calc_window(half,n,window);
    w = [w; w(end:-1:1)];
else
    % Odd length window
    half = (n+1)/2;
    w = calc_window(half,n,window);
    w = [w; w(end-1:-1:1)];
end

%---------------------------------------------------------------------
function w = calc_window(m,n,window)
%CALC_WINDOW   Calculate the generalized cosine window samples.
%   CALC_WINDOW Calculates and returns the first M points of an N point
%   generalized cosine window determined by the 'window' string.

% For the hamming and blackman windows we force rounding in order to achieve
% better numerical properties.  For example, the end points of the hamming 
% window should be exactly 0.08.

switch window
case 'hann'
    % Hann window
    %    w = 0.5 * (1 - cos(2*pi*(0:m-1)'/(n-1))); 
    a0 = 0.5;
    a1 = 0.5;
    a2 = 0;
    a3 = 0;
    a4 = 0;
case 'hamming'
    % Hamming window
    %    w = (54 - 46*cos(2*pi*(0:m-1)'/(n-1)))/100;
    a0 = 0.54;
    a1 = 0.46;
    a2 = 0;
    a3 = 0;
    a4 = 0;
case 'blackman'
    % Blackman window
    %    w = (42 - 50*cos(2*pi*(0:m-1)/(n-1)) + 8*cos(4*pi*(0:m-1)/(n-1)))'/100;
    a0 = 0.42;
    a1 = 0.5;
    a2 = 0.08;
    a3 = 0;
    a4 = 0;
case 'flattopwin'
    % Flattop window
    % Original coefficients as defined in the reference (see flattopwin.m);
    % a0 = 1;
    % a1 = 1.93;
    % a2 = 1.29;
    % a3 = 0.388;
    % a4 = 0.032;
    %
    % Scaled by (a0+a1+a2+a3+a4)
    a0 = 0.2156;
    a1 = 0.4160;
    a2 = 0.2781;
    a3 = 0.0836;
    a4 = 0.0069;
end

x = (0:m-1)'/(n-1);
w = a0 - a1*cos(2*pi*x) + a2*cos(4*pi*x) - a3*cos(6*pi*x) + a4*cos(8*pi*x);

% [EOF] gencoswin.m

%*****************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function w = hanning(varargin)
%HANNING   Hanning window.
%   HANNING(N) returns the N-point symmetric Hanning window in a column
%   vector.  Note that the first and last zero-weighted window samples
%   are not included.
%
%   HANNING(N,'symmetric') returns the same result as HANNING(N).
%
%   HANNING(N,'periodic') returns the N-point periodic Hanning window,
%   and includes the first zero-weighted window sample.
%
%   NOTE: Use the HANN function to get a Hanning window which has the 
%          first and last zero-weighted samples. 
%
%   See also BARTLETT, BLACKMAN, BOXCAR, CHEBWIN, HAMMING, HANN, KAISER
%   and TRIANG.

%   Copyright 1988-2004 The MathWorks, Inc.
%   $Revision: 1.11.4.2 $  $Date: 2004/12/26 22:16:01 $

% Check number of inputs
error(nargchk(1,2,nargin));

% Check for trivial order
[n,w,trivialwin] = check_order(varargin{1});
if trivialwin, return, end

% Select the sampling option
if nargin == 1,
   sflag = 'symmetric';
else
   sflag = lower(varargin{2});
end

% Allow partial strings for sampling options
allsflags = {'symmetric','periodic'};
sflagindex = find(strncmp(sflag, allsflags, length(sflag)));
if length(sflagindex)~=1         % catch 0 or 2 matches
   error('Sampling flag must be either ''symmetric'' or ''periodic''.');
end
sflag = allsflags{sflagindex};

% Evaluate the window
switch sflag,
case 'periodic'
   % Includes the first zero sample
   w = [0; sym_hanning(n-1)];
case 'symmetric'
   % Does not include the first and last zero sample
   w = sym_hanning(n);
end

%---------------------------------------------------------------------
function w = sym_hanning(n)
%SYM_HANNING   Symmetric Hanning window. 
%   SYM_HANNING Returns an exactly symmetric N point window by evaluating
%   the first half and then flipping the same samples over the other half.

if ~rem(n,2)
   % Even length window
   half = n/2;
   w = calc_hanning(half,n);
   w = [w; w(end:-1:1)];
else
   % Odd length window
   half = (n+1)/2;
   w = calc_hanning(half,n);
   w = [w; w(end-1:-1:1)];
end

%---------------------------------------------------------------------
function w = calc_hanning(m,n)
%CALC_HANNING   Calculates Hanning window samples.
%   CALC_HANNING Calculates and returns the first M points of an N point
%   Hanning window.

w = .5*(1 - cos(2*pi*(1:m)'/(n+1))); 

% [EOF] hanning.m

%*****************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [n_out, w, trivalwin] = check_order(n_in)
%CHECK_ORDER Checks the order passed to the window functions.
% [N,W,TRIVALWIN] = CHECK_ORDER(N_ESTIMATE) will round N_ESTIMATE to the
% nearest integer if it is not alreay an integer. In special cases (N is [],
% 0, or 1), TRIVALWIN will be set to flag that W has been modified.

%   Copyright 1988-2002 The MathWorks, Inc.
%   $Revision: 1.6 $  $Date: 2002/04/15 01:07:36 $

w = [];
trivalwin = 0;

% Special case of negative orders:
if n_in < 0,
   error('Order cannot be less than zero.');
end

% Check if order is already an integer or empty
% If not, round to nearest integer.
if isempty(n_in) | n_in == floor(n_in),
   n_out = n_in;
else
   n_out = round(n_in);
   warning('Rounding order to nearest integer.');
end

% Special cases:
if isempty(n_out) | n_out == 0,
   w = zeros(0,1);               % Empty matrix: 0-by-1
   trivalwin = 1; 
elseif n_out == 1,
   w = 1;
   trivalwin = 1;   
end
