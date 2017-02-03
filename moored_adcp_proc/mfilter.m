function [y,c]=mfilter(x,f0,f1,f2,n,m);

% MFILTER  FIR digital filter
%
% [Y,C] = MFILTER(X,F0,F1,F2,N,M);
%
% filter data vector X by N'th order FIR1 filter 
% Hamming window second order used X timeseries along 2. Dimension !!!
%
% F0 sample frequency
%
% F1 , F2 filter frequency
%
% F1=0       F2=highpass
% F1=lowpass F2=0
%
% F1 < F2 bandpass
% F1 > F2 bandstopp
%
% N Number of elements to filter, default: F0/max(F1,F2)
%
% M Increment for Subsample X after filter, default: 1
%
% C Filter-Coefficients by FIR1 with length  N+1
%
%
% see also: BFILTER, MEANIND1, MEANIND2, NOISE
%           FIR1, FILTFILT (SignalProcessingToolbox)
%
%-------------------------------------------------------------
% Example for LowPass:
%
%  x = linspace(0,50,1000);
%  y = sin(x);               % RowVector !!!
%
%  dx  = median(diff(x));    % Inkrement
%  int = 4*pi;
%
%  c   = round(int/dx);
%  z   = mfilter( y , 1 , 1/c , 0 , c , 1 );
%
%  figure,plot(x,y),hold on,plot(x,z,'r')
%

%M. Visbeck 29.07.91

%---------------------------------------------------------------------
% SubFunctions from Signal and Matlab-Toolboxes
%---------------------------------------------------------------------
%
% FILTFILT (FILT2) Zero-phase forward and reverse digital filtering.
% y = filtfilt(b,a,x)
%
% FIR1   FIR filter design using the window method.
% [b,a] = fir1(N,Wn,varargin)
%
% FIRCHK   Check if specified filter order is valid.
% [n,msg1,msg2] = firchk(n,Fend,a,exception)
%
% FIRLS Linear-phase FIR filter design using least-squares error minimization.
% [h,a]=firls(N,F,M,W,ftype);
%
% SINC Sin(pi*x)/(pi*x) function.
% y=sinc(x)
%
% HAMMING   Hamming window.
% w = hamming(varargin)
%
% GENCOSWIN   Returns one of the generalized cosine windows.
% [w,msg] = gencoswin(varargin)
%
%---------------------------------------------------------------------

if nargin < 5, n=ceil(f0/max(f1,f2)); end
if nargin < 5, m=1; end

%lp
if     f2==0, Wn = f1;      typ=''; 
%hp
elseif f1==0, Wn = f2;      typ='high'; 
%bp
elseif f2>f1, Wn = [f1 f2]; typ=''; 
%bs
elseif f1>f2, Wn = [f2 f1]; typ='stop'; 
end

Wn = 2 * Wn / f0;

b=fir1(n,Wn,typ);

si=size(x);

if any(si([1 2])) == 1
   y=filt2(b,1,x);
   y=y(1:m:end);
else
   y=x;
   for ii=1:si(1)
       y(ii,:)=filt2(b,1,x(ii,:));
   end
   y=y(:,1:m:si(2));
end

%**************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function y = filt2(b,a,x)
%FILT2 Zero-phase forward and reverse digital filtering.
%   Y = FILT2(B, A, X) filters the data in vector X with the filter described
%   by vectors A and B to create the filtered data Y.  The filter is described 
%   by the difference equation:
%
%     y(n) = b(1)*x(n) + b(2)*x(n-1) + ... + b(nb+1)*x(n-nb)
%                      - a(2)*y(n-1) - ... - a(na+1)*y(n-na)
%
%
%   After filtering in the forward direction, the filtered sequence is then 
%   reversed and run back through the filter; Y is the time reverse of the 
%   output of the second filtering operation.  The result has precisely zero 
%   phase distortion and magnitude modified by the square of the filter's 
%   magnitude response.  Care is taken to minimize startup and ending 
%   transients by matching initial conditions.
%
%   The length of the input x must be more than three times
%   the filter order, defined as max(length(b)-1,length(a)-1).
%
%   Note that FILT2 should not be used with differentiator and Hilbert FIR
%   filters, since the operation of these filters depends heavily on their
%   phase response.
%
%   See also FILTER.

%   References: 
%     [1] Sanjit K. Mitra, Digital Signal Processing, 2nd ed., McGraw-Hill, 2001
%     [2] Fredrik Gustafsson, Determining the initial states in forward-backward 
%         filtering, IEEE Transactions on Signal Processing, pp. 988--992, April 1996, 
%         Volume 44, Issue 4

%   Author(s): L. Shure, 5-17-88
%   revised by T. Krauss, 1-21-94
%   Initial Conditions: Fredrik Gustafsson
%   Copyright 1988-2004 The MathWorks, Inc.
%   $Revision: 1.7.4.2 $  $Date: 2004/12/26 22:15:49 $

    error(nargchk(3,3,nargin))
    if (isempty(b) || isempty(a) || isempty(x))
        y = [];
        return
    end

    [m,n] = size(x);
    if (n>1) && (m>1)
        y = x;
        for i=1:n  % loop over columns
            y(:,i) = filt2(b,a,x(:,i));
        end
        return
        % error('Only works for vector input.')
    end
    if m==1
        x = x(:);   % convert row to column
    end
    len = size(x,1);   % length of input
    b = b(:).';
    a = a(:).';
    nb = length(b);
    na = length(a);
    nfilt = max(nb,na);

    nfact = 3*(nfilt-1);  % length of edge transients

    if (len<=nfact),    % input data too short!
        error('Data must have length more than 3 times filter order.');
    end

% set up filter's initial conditions to remove dc offset problems at the 
% beginning and end of the sequence
    if nb < nfilt, b(nfilt)=0; end   % zero-pad if necessary
    if na < nfilt, a(nfilt)=0; end
% use sparse matrix to solve system of linear equations for initial conditions
% zi are the steady-state states of the filter b(z)/a(z) in the state-space 
% implementation of the 'filter' command.
    rows = [1:nfilt-1  2:nfilt-1  1:nfilt-2];
    cols = [ones(1,nfilt-1) 2:nfilt-1  2:nfilt-1];
    data = [1+a(2) a(3:nfilt) ones(1,nfilt-2)  -ones(1,nfilt-2)];
    sp = sparse(rows,cols,data);
    zi = sp \ ( b(2:nfilt).' - a(2:nfilt).'*b(1) );
% non-sparse:
% zi = ( eye(nfilt-1) - [-a(2:nfilt).' [eye(nfilt-2); zeros(1,nfilt-2)]] ) \ ...
%      ( b(2:nfilt).' - a(2:nfilt).'*b(1) );

% Extrapolate beginning and end of data sequence using a "reflection
% method".  Slopes of original and extrapolated sequences match at
% the end points.
% This reduces end effects.
    y = [2*x(1)-x((nfact+1):-1:2);x;2*x(len)-x((len-1):-1:len-nfact)];

% filter, reverse data, filter again, and reverse data again
    y = filter(b,a,y,zi*y(1));
    y = y(length(y):-1:1);
    y = filter(b,a,y,zi*y(1));
    y = y(length(y):-1:1);

% remove extrapolated pieces of y
    y([1:nfact len+nfact+(1:nfact)]) = [];

    if m == 1
        y = y.';   % convert back to row if necessary
    end


%**************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [b,a] = fir1(N,Wn,varargin)
%FIR1   FIR filter design using the window method.
%   B = FIR1(N,Wn) designs an N'th order lowpass FIR digital filter
%   and returns the filter coefficients in length N+1 vector B.
%   The cut-off frequency Wn must be between 0 < Wn < 1.0, with 1.0 
%   corresponding to half the sample rate.  The filter B is real and
%   has linear phase.  The normalized gain of the filter at Wn is
%   -6 dB.
%
%   B = FIR1(N,Wn,'high') designs an N'th order highpass filter.
%   You can also use B = FIR1(N,Wn,'low') to design a lowpass filter.
%
%   If Wn is a two-element vector, Wn = [W1 W2], FIR1 returns an
%   order N bandpass filter with passband  W1 < W < W2. You can
%   also specify B = FIR1(N,Wn,'bandpass').  If Wn = [W1 W2], 
%   B = FIR1(N,Wn,'stop') will design a bandstop filter.
%
%   If Wn is a multi-element vector, 
%          Wn = [W1 W2 W3 W4 W5 ... WN],
%   FIR1 returns an order N multiband filter with bands
%    0 < W < W1, W1 < W < W2, ..., WN < W < 1.
%   B = FIR1(N,Wn,'DC-1') makes the first band a passband.
%   B = FIR1(N,Wn,'DC-0') makes the first band a stopband.
%
%   B = FIR1(N,Wn,WIN) designs an N-th order FIR filter using 
%   the N+1 length vector WIN to window the impulse response.
%   If empty or omitted, FIR1 uses a Hamming window of length N+1.
%   For a complete list of available windows, see the help for the
%   WINDOW function. KAISER and CHEBWIN can be specified with an 
%   optional trailing argument.  For example, B = FIR1(N,Wn,kaiser(N+1,4)) 
%   uses a Kaiser window with beta=4. B = FIR1(N,Wn,'high',chebwin(N+1,R)) 
%   uses a Chebyshev window with R decibels of relative sidelobe 
%   attenuation.
%
%   For filters with a gain other than zero at Fs/2, e.g., highpass
%   and bandstop filters, N must be even.  Otherwise, N will be
%   incremented by one.  In this case the window length should be
%   specified as N+2.
%   
%   By default, the filter is scaled so the center of the first pass band 
%   has magnitude exactly one after windowing. Use a trailing 'noscale' 
%   argument to prevent this scaling, e.g. B = FIR1(N,Wn,'noscale'), 
%   B = FIR1(N,Wn,'high','noscale'), B = FIR1(N,Wn,wind,'noscale').  You
%   can also specify the scaling explicitly, e.g. FIR1(N,Wn,'scale'), etc.
%
%   See also KAISERORD, FIRCLS1, FIR2, FIRLS, FIRCLS, CFIRPM,
%            FIRPM, FREQZ, FILTER, WINDOW.

%   FIR1 is an M-file implementation of program 5.2 in the IEEE
%   Programs for Digital Signal Processing tape. 

%   Author(s): L. Shure
%              L. Shure, 4-5-90, revised
%              T. Krauss, 3-5-96, revised
%   Copyright 1988-2004 The MathWorks, Inc.
%   $Revision: 1.15.4.2 $  $Date: 2004/04/13 00:17:49 $

%   Reference(s):
%     [1] "Programs for Digital Signal Processing", IEEE Press
%         John Wiley & Sons, 1979, pg. 5.2-1.

error(nargchk(2,5,nargin));

% Parse optional input arguments
[Ftype,Wind,SCALING,msg] = parseoptargs(Wn,varargin{:});
error(msg);

% Compute the frequency vector
[nbands,ff,Ftype,msg] = desiredfreq(Wn,Ftype);
error(msg);

% Compute the magnitude vector
[aa,First_Band] = desiredmag(Ftype,nbands);

% Check for appropriate filter order, increase when necessary
[N,msg1,msg2] = firchk(N,ff(end),aa);
error(msg1);
warning(msg2);

% Work with filter length (= order + 1)
L = N + 1;

% Check for valid window, or assign default if empty
[Wind,msg] = chkwindow(Wind,L);
error(msg);

% Compute unwindowed impulse response
hh = firls(L-1,ff,aa);

% Window impulse response to get the filter
b = hh.*Wind(:)'; 
a = 1;

if SCALING,   
    % Scale so that passband is approx 1
    b = scalefilter(b,First_Band,ff,L);   
end

%-------------------------------------------------------------------------
function [Ftype,Wind,SCALING,msg] = parseoptargs(Wn,varargin)
%PARSEOPTARGS   Parse optional input arguments.

% Up to 3 optional input arguments, always in this order:
%   1 - Filter type flag, can be 'low','high','bandpass','stop','DC-0','DC-1'
%   2 - Window vector
%   3 - 'noscale' flag

% Initialize output args.
msg = '';
SCALING = [];

[Ftype,Wind,Scale] = assignoptargs(Wn,varargin{:});

[Ftype,Wind,Scale,msg] = validateargs(Wn,Ftype,Wind,Scale);
if ~isempty(msg),
    return
end

switch lower(Scale),
case 'noscale';
    SCALING = 0;
case 'scale';
    SCALING = 1;
end

%--------------------------------------------------------------------------
function [Ftype,Wind,Scale] = assignoptargs(Wn,varargin)
%ASSIGNOPTARGS  Assign optional input arguments to the appropriate variables.

% default optional parameter values:
Wind = [];
Scale = 'scale';
Ftype = defaultftype(Wn);


switch length(varargin)
case 1
    if ischar(varargin{1}) && (length(varargin{1})>0),
        s = upper(varargin{1});
        switch upper(s)
        case {'SCALE','NOSCALE'}
            Scale = s;
        otherwise
            Ftype = s;
        end
    else
        Wind = varargin{1};
    end
case 2
    if ischar(varargin{1})
        Ftype = varargin{1};
    else
        Wind = varargin{1};
    end
    if ischar(varargin{2})
        Scale = varargin{2};
    else
        Wind = varargin{2};
    end
case 3
    Ftype = varargin{1};
    Wind = varargin{2};
    Scale = varargin{3};
end

%--------------------------------------------------------------------------
function [Ftype,Wind,Scale,msg] = validateargs(Wn,Ftype,Wind,Scale)
%VALIDATEARGS  Test if arguments are valid.

msg = '';

% Assign a default Ftype when an empty is given. Backwards compatibility
if isempty(Ftype),
    Ftype = defaultftype(Wn);
end


Ftypeopts = {'LOW','HIGH','BANDPASS','STOP','DC-0','DC-1'};
Scaleopts = {'NOSCALE','SCALE'};


indx = strmatch(upper(Ftype),Ftypeopts);
if isempty(indx),
    msg = 'Unrecognized or ambiguous filter type specified.';
    return
else
    Ftype = Ftypeopts{indx};
end

scaleindx = strmatch(upper(Scale),Scaleopts);
if isempty(scaleindx),
    msg = 'Scaling option must be ''noscale'' or ''scale''.';
    return
else
    Scale = Scaleopts{scaleindx};
end

if ~any(size(Wind) <= 1),
    msg = 'The window specified must be a vector';
    return
else
    Wind = Wind(:).'; % Make it a row vector
end

%--------------------------------------------------------------------------
function [nbands,ff,Ftype,msg] = desiredfreq(Wn,Ftype)
%DESIREDFREQ  Compute the vector of frequencies to pass to FIRLS.
%
%   Inputs:
%           Wn    - vector of cutoff frequencies.
%           Ftype - string with desired response ('low','high',...)
%
%   Outputs:
%           nbands - number of frequency bands.
%           ff     - vector of frequencies to pass to FIRLS.
%           Ftype  - converted filter type (if it's necessary to convert)

% Initialize output args.
nbands = [];
ff     = [];
msg    = '';


if  any( Wn<0 | Wn>1 ),
   msg = 'Frequencies must fall in range between 0 and 1.';
   return
end
if  any(diff(Wn)<0),
   msg = 'Frequencies must be increasing';
   return
end

Wn = Wn(:)';

nbands = length(Wn) + 1;

if (nbands > 2) && strcmp(lower(Ftype),'bandpass'),
    Ftype = 'DC-0';  % make sure default 3 band filter is bandpass
end

ff = [0,Wn(1:nbands-1); Wn(1:nbands-1),1];

ff = ff(:);

%-------------------------------------------------------------------------
function [aa,First_Band] = desiredmag(Ftype,nbands)
%DESIREDMAG  Compute the magnitude vector to pass to FIRLS.

First_Band = isempty(findstr('DC-0',Ftype)) && isempty(findstr('HIGH',Ftype));
mags = rem( First_Band + (0:nbands-1), 2);
aa = [mags(:)'; mags(:)'];

aa = aa(:);
%--------------------------------------------------------------------------
function [Wind,msg] = chkwindow(Wind,L)
%CHKWINDOW   Check if specified window is valid, assign default if empty.

msg = '';

if isempty(Wind),
   % Replace the following with the default window of your choice.
   Wind = hamming(L);
end
    
if length(Wind) ~= L
   msg = 'The window length must be the same as the filter length.';
end
%
% to use Kaiser window, beta must be supplied
% att = 60; % dB of attenuation desired in sidelobe
% beta = 0.1102*(att-8.7);
% wind = kaiser(L,beta);

%---------------------------------------------------------------------------
function b = scalefilter(b,First_Band,ff,L)
%SCALEFILTER   Scale fitler to have passband approx. equal to one.

if First_Band
    b = b / sum(b);  % unity gain at DC
else
    if ff(4)==1
        % unity gain at Fs/2
        f0 = 1;
    else
        % unity gain at center of first passband
        f0 = mean(ff(3:4));
    end
    b = b / abs( exp(-j*2*pi*(0:L-1)*(f0/2))*(b.') );
end

%----------------------------------------------------------------------------
function Ftype = defaultftype(Wn)
%DEFAULTFTYPE  Assign default filter type depending on number of bands.

if length(Wn) == 1,
    Ftype = 'low';
elseif length(Wn) == 2,
    Ftype = 'bandpass';
elseif length(Wn) >= 3,
    Ftype = 'dc-0';
end

%**************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [n,msg1,msg2] = firchk(n,Fend,a,exception)
%FIRCHK   Check if specified filter order is valid.
%   FIRCHK(N,Fend,A) checks if the specified order N is valid given the
%   final frequency point Fend and the desired magnitude response vector A.
%   Type 2 linear phase FIR filters (symmetric, odd order) must have a
%   desired magnitude response vector that ends in zero if Fend = 1.  This
%   is because type 2 filters necessarily have a zero at w = pi.
%
%   If the order is not valid, a warning is given and the order
%   of the filter is incremented by one.
%
%   If A is a scalar (as when called from fircls1), A = 0 is
%   interpreted as lowpass and A = 1 is interpreted as highpass.
%
%   FIRCHK(N,Fend,A,EXCEPTION) will not warn or increase the order
%   if EXCEPTION = 1.  Examples of EXCEPTIONS are type 4 filters
%   (such as differentiators or hilbert transformers) or non-linear
%   phase filters (such as minimum and maximum phase filters).

%   Author : R. Losada
%   Copyright 1988-2004 The MathWorks, Inc.
%   $Revision: 1.7.4.3 $  $Date: 2004/04/13 00:18:44 $

nargchk(3,4,nargin);

if nargin == 3,
    exception = false;
end

msg1 = '';
msg2 = '';
oddord = false; % Flag, initially we assume even order

if isempty(n) || length(n) > 1 || ~isnumeric(n) || ~isreal(n) || n~=round(n) || n<=0,
    msg1 = 'Filter order must be a real, positive integer.';
    return
end

if rem(n,2) == 1,
    oddord = true; % Overwrite flag
end
 
if (a(end) ~= 0) && Fend == 1 && oddord && ~exception,
    str = ['Odd order symmetric FIR filters must have a gain of zero \n'...
     'at the Nyquist frequency. The order is being increased by one.'];
    msg2 = sprintf(str);
    n = n+1;
end
    
%**************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [h,a]=firls(N,F,M,W,ftype);
% FIRLS Linear-phase FIR filter design using least-squares error minimization.
%   B=FIRLS(N,F,A) returns a length N+1 linear phase (real, symmetric
%   coefficients) FIR filter which has the best approximation to the
%   desired frequency response described by F and A in the least squares
%   sense. F is a vector of frequency band edges in pairs, in ascending
%   order between 0 and 1. 1 corresponds to the Nyquist frequency or half
%   the sampling frequency. A is a real vector the same size as F
%   which specifies the desired amplitude of the frequency response of the
%   resultant filter B. The desired response is the line connecting the
%   points (F(k),A(k)) and (F(k+1),A(k+1)) for odd k; FIRLS treats the
%   bands between F(k+1) and F(k+2) for odd k as "transition bands" or
%   "don't care" regions. Thus the desired amplitude is piecewise linear
%   with transition bands.  The integrated squared error is minimized.
%
%   For filters with a gain other than zero at Fs/2, e.g., highpass
%   and bandstop filters, N must be even.  Otherwise, N will be
%   incremented by one. Alternatively, you can use a trailing 'h' flag to
%   design a type 4 linear phase filter and avoid incrementing N.
%
%   B=FIRLS(N,F,A,W) uses the weights in W to weight the error. W has one
%   entry per band (so it is half the length of F and A) which tells
%   FIRLS how much emphasis to put on minimizing the integral squared error
%   in each band relative to the other bands.
%
%   B=FIRLS(N,F,A,'Hilbert') and B=FIRLS(N,F,A,W,'Hilbert') design filters
%   that have odd symmetry, that is, B(k) = -B(N+2-k) for k = 1, ..., N+1.
%   A special case is a Hilbert transformer which has an approx. amplitude
%   of 1 across the entire band, e.g. B=FIRLS(30,[.1 .9],[1 1],'Hilbert').
%
%   B=FIRLS(N,F,A,'differentiator') and B=FIRLS(N,F,A,W,'differentiator')
%   also design filters with odd symmetry, but with a special weighting
%   scheme for non-zero amplitude bands. The weight is assumed to be equal
%   to the inverse of frequency, squared, times the weight W. Thus the
%   filter has a much better fit at low frequency than at high frequency.
%   This designs FIR differentiators.
%
%   % Example of a length 31 lowpass filter:
%   	h=firls(30,[0 .1 .2 .5]*2,[1 1 0 0]);
%
%   % Example of a low-pass differentiator:
%   	h=firls(44,[0 .3 .4 1],[0 .2 0 0],'differentiator');
%
%   % Example of a type 4 highpass filter:
%       h=firls(25,[0 .4 .5 1],[0 0 1 1],'h');
%
%   See also FIRPM, FIR1, FIR2, FREQZ and FILTER.

%       Author(s): T. Krauss
%   History: 10-18-91, original version
%            3-30-93, updated
%            9-1-95, optimize adjacent band case
%   Copyright 1988-2004 The MathWorks, Inc.
%   $Revision: 1.11.4.2 $  $Date: 2004/04/13 00:17:54 $

% check number of arguments, set up defaults.
nargchk(3,5,nargin);

if (max(F)>1) || (min(F)<0)
    error('Frequencies in F must be in range [0,1].')
end
if (rem(length(F),2)~=0)
    error('F must have even length.');
end
if (length(F) ~= length(M))
    error('F and A must be equal lengths.');
end
if (nargin==3),
    W = ones(length(F)/2,1);
    ftype = '';
end
if (nargin==4),
    if isstr(W),
        ftype = W; W = ones(length(F)/2,1);
    else
        ftype = '';
    end
end
if (nargin==5),
    if isempty(W),
        W = ones(length(F)/2,1);
    end
end
if isempty(ftype)
    ftype = 0;  differ = 0;
else
    ftype = lower(ftype);
    if strcmpi(ftype,'h') || strcmpi(ftype,'hilbert')
        ftype = 1;  differ = 0;
    elseif strcmpi(ftype,'d') || strcmpi(ftype,'differentiator')
        ftype = 1;  differ = 1;
    else
        error('Requires symmetry to be ''Hilbert'' or ''differentiator''.')
    end
end

% Check for valid filter length
[N,msg1,msg2] = firchk(N,F(end),M,ftype);
error(msg1);

if ~isempty(msg2),
    msg2 = sprintf([msg2,'\r',...
        '\nAlternatively, you can pass a trailing ''h'' argument,\r',...
        'as in firls(N,F,A,W,''h''), to design a type 4 linear phase filter.']);
end
warning(msg2);


N = N+1;                   % filter length
F=F(:)/2;  M=M(:);  W=sqrt(W(:));  % make these guys columns
dF = diff(F);

if (length(F) ~= length(W)*2)
    error('There should be one weight per band.');
end;
if any(dF<0),
    error('Frequencies in F must be nondecreasing.')
end
if all(dF(2:2:length(dF)-1)==0)
    fullband = 1;
else
    fullband = 0;
end
if all((W-W(1))==0)
    constant_weights = 1;
else
    constant_weights = 0;
end

L=(N-1)/2;

Nodd = rem(N,2);

if (ftype == 0),  % Type I and Type II linear phase FIR
    % basis vectors are cos(2*pi*m*f) (see m below)
    if ~Nodd
        m=(0:L)+.5;   % type II
    else
        m=(0:L);      % type I
    end
    k=m';
    need_matrix = (~fullband) || (~constant_weights);
    if need_matrix
        I1=k(:,ones(size(m)))+m(ones(size(k)),:);    % entries are m + k
        I2=k(:,ones(size(m)))-m(ones(size(k)),:);    % entries are m - k
        G=zeros(size(I1));
    end

    if Nodd
        k=k(2:length(k));
        b0=0;       %  first entry must be handled separately (where k(1)=0)
    end;
    b=zeros(size(k));
    for s=1:2:length(F),
        m=(M(s+1)-M(s))/(F(s+1)-F(s));    %  slope
        b1=M(s)-m*F(s);                   %  y-intercept
        if Nodd
            b0 = b0 + (b1*(F(s+1)-F(s)) + m/2*(F(s+1)*F(s+1)-F(s)*F(s)))...
                * abs(W((s+1)/2)^2) ;
        end
        b = b+(m/(4*pi*pi)*(cos(2*pi*k*F(s+1))-cos(2*pi*k*F(s)))./(k.*k))...
            * abs(W((s+1)/2)^2);
        b = b + (F(s+1)*(m*F(s+1)+b1)*sinc(2*k*F(s+1)) ...
            - F(s)*(m*F(s)+b1)*sinc(2*k*F(s))) ...
            * abs(W((s+1)/2)^2);
        if need_matrix
            G = G + (.5*F(s+1)*(sinc(2*I1*F(s+1))+sinc(2*I2*F(s+1))) ...
                - .5*F(s)*(sinc(2*I1*F(s))+sinc(2*I2*F(s))) ) ...
                * abs(W((s+1)/2)^2);
        end
    end;
    if Nodd
        b=[b0; b];
    end;

    if need_matrix
        a=G\b;
    else
        a=(W(1)^2)*4*b;
        if Nodd
            a(1) = a(1)/2;
        end
    end
    if Nodd
        h=[a(L+1:-1:2)/2; a(1); a(2:L+1)/2].';
    else
        h=.5*[flipud(a); a].';
    end;
elseif (ftype == 1),  % Type III and Type IV linear phase FIR
    %  basis vectors are sin(2*pi*m*f) (see m below)
    if (differ),      % weight non-zero bands with 1/f^2
        do_weight = ( abs(M(1:2:length(M))) +  abs(M(2:2:length(M))) ) > 0;
    else
        do_weight = zeros(size(F));
    end

    if Nodd
        m=(1:L);      % type III
    else
        m=(0:L)+.5;   % type IV
    end;
    k=m';
    b=zeros(size(k));

    need_matrix = (~fullband) || (any(do_weight)) || (~constant_weights);
    if need_matrix
        I1=k(:,ones(size(m)))+m(ones(size(k)),:);    % entries are m + k
        I2=k(:,ones(size(m)))-m(ones(size(k)),:);    % entries are m - k
        G=zeros(size(I1));
    end

    i = sqrt(-1);
    for s=1:2:length(F),
        if (do_weight((s+1)/2)),      % weight bands with 1/f^2
            if F(s) == 0, F(s) = 1e-5; end     % avoid singularities
            m=(M(s+1)-M(s))/(F(s+1)-F(s));
            b1=M(s)-m*F(s);
            snint1 = sineint(2*pi*k*F(s+1)) - sineint(2*pi*k*F(s));
            %snint1 = (-1/2/i)*(expint(i*2*pi*k*F(s+1)) ...
            %    -expint(-i*2*pi*k*F(s+1)) -expint(i*2*pi*k*F(s)) ...
            %    +expint(-i*2*pi*k*F(s)) );
            % csint1 = cosint(2*pi*k*F(s+1)) - cosint(2*pi*k*F(s)) ;
            csint1 = (-1/2)*(expint(i*2*pi*k*F(s+1))+expint(-i*2*pi*k*F(s+1))...
                -expint(i*2*pi*k*F(s))  -expint(-i*2*pi*k*F(s)) );
            b=b + ( m*snint1 ...
                + b1*2*pi*k.*( -sinc(2*k*F(s+1)) + sinc(2*k*F(s)) + csint1 ))...
                * abs(W((s+1)/2)^2);
            snint1 = sineint(2*pi*F(s+1)*(-I2));
            snint2 = sineint(2*pi*F(s+1)*I1);
            snint3 = sineint(2*pi*F(s)*(-I2));
            snint4 = sineint(2*pi*F(s)*I1);
            G = G - ( ( -1/2*( cos(2*pi*F(s+1)*(-I2))/F(s+1)  ...
                - 2*snint1*pi.*I2 ...
                - cos(2*pi*F(s+1)*I1)/F(s+1) ...
                - 2*snint2*pi.*I1 )) ...
                - ( -1/2*( cos(2*pi*F(s)*(-I2))/F(s)  ...
                - 2*snint3*pi.*I2 ...
                - cos(2*pi*F(s)*I1)/F(s) ...
                - 2*snint4*pi.*I1) ) ) ...
                * abs(W((s+1)/2)^2);
        else      % use usual weights
            m=(M(s+1)-M(s))/(F(s+1)-F(s));
            b1=M(s)-m*F(s);
            b=b+(m/(4*pi*pi)*(sin(2*pi*k*F(s+1))-sin(2*pi*k*F(s)))./(k.*k))...
                * abs(W((s+1)/2)^2) ;
            b = b + (((m*F(s)+b1)*cos(2*pi*k*F(s)) - ...
                (m*F(s+1)+b1)*cos(2*pi*k*F(s+1)))./(2*pi*k)) ...
                * abs(W((s+1)/2)^2) ;
            if need_matrix
                G = G + (.5*F(s+1)*(sinc(2*I1*F(s+1))-sinc(2*I2*F(s+1))) ...
                    - .5*F(s)*(sinc(2*I1*F(s))-sinc(2*I2*F(s)))) * ...
                    abs(W((s+1)/2)^2);
            end
        end;
    end

    if need_matrix
        a=G\b;
    else
        a=-4*b*(W(1)^2);
    end
    if Nodd
        h=.5*[flipud(a); 0; -a].';
    else
        h=.5*[flipud(a); -a].';
    end
    if differ, h=-h; end
end

if nargout > 1
    a = 1;
end

%----------------------------------------------------------------------------
function y = sineint(x)
% SINEINT (a.k.a. SININT)   Numerical Sine Integral
%   Used by FIRLS in the Signal Processing Toolbox.
%   Untested for complex or imaginary inputs.
%
%   See also SININT in the Symbolic Toolbox.

%   Was Revision: 1.5, Date: 1996/03/15 20:55:51

i1 = find(real(x)<0);   % this equation is not valid if x is in the
% left-hand plane of the complex plane.
% use relation Si(-z) = -Si(z) in this case (Eq 5.2.19, Abramowitz
%  & Stegun).
x(i1) = -x(i1);
y = zeros(size(x));
ind = find(x);
% equation 5.2.21 Abramowitz & Stegun
%  y(ind) = (1/(2*i))*(expint(i*x(ind)) - expint(-i*x(ind))) + pi/2;
y(ind) = imag(expint(i*x(ind))) + pi/2;
y(i1) = -y(i1);

%**************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function y=sinc(x)
%SINC Sin(pi*x)/(pi*x) function.
%   SINC(X) returns a matrix whose elements are the sinc of the elements 
%   of X, i.e.
%        y = sin(pi*x)/(pi*x)    if x ~= 0
%          = 1                   if x == 0
%   where x is an element of the input matrix and y is the resultant
%   output element.
%
%   See also SQUARE, SIN, COS, CHIRP, DIRIC, GAUSPULS, PULSTRAN, RECTPULS,
%   and TRIPULS.

%   Author(s): T. Krauss, 1-14-93
%   Copyright 1988-2002 The MathWorks, Inc.
%       $Revision: 1.7 $  $Date: 2002/04/15 01:13:58 $

y=ones(size(x));
i=find(x);
y(i)=sin(pi*x(i))./(pi*x(i));

%**************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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

%**************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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

%----------------------------------------------------------------------------

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
