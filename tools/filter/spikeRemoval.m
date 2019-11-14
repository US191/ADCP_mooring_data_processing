function [x,spikes] = spikeRemoval(w,varargin)
%   DESPIKING DISCRETE-TIME SIGNAL USING HISTOGRAM METHOD (a.k.a. DELETE OUTLIERS)
%
%   Time series may contain undesired transients and spikes. This function
%   replaces sudden spikes (outliers) exceeding the threshold value by
%   either:
%
%   (1) interpolating among previous and subsequent data points, or
%   (2) replacing them with NaN (delete outliers)
%
%   User needs to choose one of the above options; default is interpolation.
%
%   "w" is the input array and "x" is the spike removed output array. In
%   addition, output argument "spikes" is a structure array, which returns
%   the indices of spikes, corresponding values, and replacement values. If
%   npass parameter is selected as 2 (it means the code will run twice),
%   then "spikes" will be the structure array of results corresponding to
%   each run.
%
%   The threshold is defined as mean +/- a number of standard deviations of
%   windowed data centered at spike locations.
%
%   USAGE:
%   [x,spikes] = spikeRemoval(w) 
%
%   or
%   
%   [x,spikes] = spikeRemoval(w,prop_name,prop_val)
%
%   STATIC INPUT:
%            w = vector of time series data (1xn or nx1)
%
%   VALID PROP_NAME / PROP_VAL PAIRS:
%   -----------------------------------------
%   'wnsz'   --> (1x1)-[numeric]-[default:25]
%   'nstd'   --> (1x2)-[numeric]-[default: 3]
%   'nbins'  --> (1x2)-[numeric]-[default: 4]
%   'ndata'  --> (1x1)-[numeric]-[default: 5]
%   'npass'  --> (1x1)-[numeric]-[default: 2]
%   'method' --> [text]-[default: 'Linear']
%   'debug'  --> [text]-[default: 'True']
%
%  NOTES:
%         x = spike removed time series
%
%    spikes = 3xn structure array, first column contains indexes of spikes,
%             second column contains values of spikes and third column
%             contains replacement values
%
%      wnsz = window size to compute mean and standard deviation for
%             thresholding
%
%      nstd = number of standard deviations above and below mean
%
%     nbins = number of bins used in histogram method
%
%     ndata = number of neighbouring data points for spike replacement
%
%     npass = number of passes for spike detection and replacement (e.g.,
%             1); Suggest npass = 2 or 3 for noisy data with back-to-back
%             spikes
%
%    method = interpolation method (linear, spline, cubic, nearest). To
%             delete outliers chose method as 'delete' 
%             (e.g., 'method','delete').
%
%     debug = 'True' for print debug messages and plot results or '' for no
%             debug messages or plots
%
%   WARNING: This function is for detecting sudden spikes, it may not be
%            effective detecting long duration spike-like pulses.
%
%   EXAMPLE 1 (using defaults values):
%
%    w = rand(500,1)*10; w(ceil(rand(1,20)*500))=rand(1,20)*100;
%    [x,spikes] = spikeRemoval(w);
%
%   EXAMPLE 2 (using user defined values for interpolation):
%
%    w = rand(500,1)*10; w(ceil(rand(1,20)*500))=rand(1,20)*100;
%    [x,spikes] = spikeRemoval(w,'wnsz',20,'nstd',2,'npass',1);
%
%   EXAMPLE 3 (using user defined values for deleting outliers):
%
%    w = rand(500,1)*10; w(ceil(rand(1,20)*500))=rand(1,20)*100;
%    [x,spikes] = spikeRemoval(w,'nstd',2,'npass',1,'method','delete');
%
%   REQUIREMENTS:
%   spikeRemoval function doesn't require any MatLAB toolbox.
%
%   ACKNOWLEDGEMENTS:
%   I would like thank Dr. Ayal Anis of Texas AM and Jeanne Jones of USGS
%   for their valuable comments and suggestions, which helped improving
%   this code significantly. Also, thanks to Jan Glscher of University of
%   Hamburg for making his nanmean and nanstd codes available. Finally,
%   thanks to James for making his cent_diff_n function available at MatLAB
%   FEX.
%
%   REFERENCES:
%
%   Solomon, O. M., D. R. Larson, and N. G. Paulter (2001). Comparison of
%   some algorithms to estimate the low and high state level of pulses,
%   IEEE Instrumentation and Measurement Technology Conference, Vol. 1,
%   Budapest, Hungary, 21?23 May 2001, 96?101.
%
%   cent_diff_n function is available at
%   https://www.mathworks.com/matlabcentral/fileexchange/36123-n-point-central-differencing
%   
%   nanmean and nanstd functions are available at
%   https://www.mathworks.com/matlabcentral/fileexchange/6837-nan-suite
%
%   THIS SOFTWARE IS PROVIDED "AS IS" AND ANY EXPRESS OR IMPLIED
%   WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
%   MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN
%   NO EVENT SHALL THE COPYRIGHT OWNER BE LIABLE FOR ANY DIRECT, INDIRECT,
%   INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
%   BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
%   OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
%   ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR
%   TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
%   USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
%   DAMAGE.
%
%   Written by Dr. Erol Kalkan, P.E. (kalkan76@gmail.com)
%   Revision: 2.1.0
%   Date: 2019/03/20 08:58:00 $

%% DEFAULT PROPERTIES
if (nargin >= 1)
    wnsz = 25;
    nstd = 3;
    ndata = 5;
    nbins = 4;
    method = 'Linear';
    debug = 'True';
    npass = 2;
else
    error('spikeRemoval: First argument must be a waveform')
end

%% USER-DEFINED PROPERTIES
if (nargin > 1)
    v = varargin;
    nv = nargin-1;
    if ~rem(nv,2) == 0
        error(['spikeRemoval: Arguments after wave must appear in ',...
            'property name/val pairs'])
    end
    for n = 1:2:nv-1
        name = lower(v{n});
        val = v{n+1};
        switch name
            case 'nstd'
                nstd = val;
            case 'wnsz'
                wnsz = val;
            case 'ndata'
                ndata = val;
            case 'method'
                method = val;
            case 'debug'
                debug = val;
            case 'nbins'
                nbins = val;
            case 'npass'
                npass = val;
            otherwise
                error('spikeRemoval: Property name not recognized')
        end
    end
end

% Initialize
nspk  = 0; % number of spikes (outliers)
x = w;

for i = 1:npass
    [x,spikes(i).results,nspk] = run(x,nspk,nstd,wnsz,ndata,method,debug,nbins);
end
if debug
    fprintf('Number of spikes found and replaced = %4d\n',nspk);
    figure(i+1);
    subplot(311); plot(w,'r'); hold on; title('Original'); lim = ylim;
    subplot(312); plot(x); title('SpikeRemoval'); ylim(lim);
    subplot(313); plot(medfilt1_perso(w,nstd),'g');
    title('MatLAB 1-D Median Filter'); ylim(lim);
end
end

function [x,spikes,nspk] = run(x,nspk,nstd,wnsz,ndata,method,debug,nbins)
% Enforce input as row vector
if ~isrow(x)
    x = x';
end

% Compute first derivative using central differentiation with order 3
xdot = cent_diff_n(x,1,3);

% Define spike locations via histogram method
R = statelevel(abs(xdot),nbins);
locs = find(abs(xdot) > R(1));

if debug
    figure(1); plot(xdot); hold on; plot(locs,xdot(locs),'ro');
    title('Spikes to be checked on first derivative of time series');
end

for i = 1:length(locs)
    [x,spike,flag] = spikeFix(x,locs(i),nstd,wnsz,ndata,method,debug);
    nspk = nspk + flag;
    if nspk ~= 0 && flag == 1
        spikes(nspk,:) = spike;
    end
end
if ~exist('spikes','var')
    spikes = -1;
end
end

function [x,spike,flag] = spikeFix(x,idx,nstd,wnsz,ndata,method,debug)
%   Input:
%         x = vector of time series (1xn)
%
%       idx = spike location
%
%   Output:
%         x = vector of spiked removed time series
%
%      flag = 1 or 0; 1 means spike is found and fixed, 0 means no spike
%             detected

% If spike is at the beginning, replace it with the mean value of next
% ndata points
if idx <= 2*ndata
    w = x(idx);
    if method == 'delete'
        x(idx) = NaN;
    else
        x(idx) = nanmean(x(idx+1:idx+ndata));
    end
    flag = 1;
    spike = [idx,w,x(idx)];
    if debug
        fprintf('Peak of %5.2f @%5.0f is replaced by %5.2f\n',w,idx,x(idx));
    end
    % If spike is at the end, replace it with the mean value of previous
    % ndata points
elseif idx >= length(x)-2*ndata
    w = x(idx);
    if method == 'delete'
        x(idx) = NaN;
    else
        x(idx) = nanmean(x(idx-ndata:idx-1));
    end
    flag = 1;
    spike = [idx,w,x(idx)];
    if debug
        fprintf('Peak of %5.2f @%5.0f is replaced by %5.2f\n',w,idx,x(idx));
    end
else
    % Set spike search region. Note: arr needs to be updated if n-point
    % central difference method is used for differentiation, for instance
    % if 5-point central difference method is used then arr =
    % (idx-2:idx+1); the current implementation uses a 3-point central
    % difference
    arr = (idx-1:idx+1);
    [~, idx] = max(abs(x(arr)));
    i = arr(idx);
    
    % Overwrite window size if spike occurs within first or last data
    % window
    if i <= floor(wnsz/2)
        wnsz = (i-1)*2;
    elseif i >= length(x) - floor(wnsz/2)
        wnsz = length(x) - i;
    end
    
    % Overwrite ndata if it is bigger than half of window size
    if ndata > floor(wnsz/2)
        ndata = floor(wnsz/2)-1;
    end
    
    % Compute mean and standard deviation of data window centered at spike
    % by excluding spike
    arr1 = [x((i - floor(wnsz/2)):(i-1)), x((i+1):(i + floor(wnsz/2)))];
    mu_x = nanmean(arr1); std_x = nanstd(arr1);
    
    if debug
        fprintf('wnsz =%5.0f, ndata =%5.0f, i = %5.0f, x(i) = %5.2f, mu_x = %5.2f, std_x = %5.2f\n', ...
            wnsz,ndata,i,x(i),mu_x,std_x);
    end
    
    % Check threshold exceedance
    if x(i) >= mu_x + std_x*nstd || x(i) <= mu_x - std_x*nstd
        w = x(i);
        
        if method == 'delete'
            x(i) = NaN;
        else
            % Replace spike region using interpolation of neighbouring data
            % points
            x(i-1:i+1) = interp1([i-ndata-1:i-2,i+2:i+ndata+1],...
                [x(i-ndata-1:i-2),x(i+2:i+ndata+1)],i-1:i+1,method);
        end
        flag = 1;
        spike = [i,w,x(i)];
        if debug
            fprintf('Peak of %5.2f @%5.0f is replaced by %5.2f\n',w,i,x(i));
        end
    else
        flag = 0;
        spike = -1;
    end
end
end

function [levels, histogram, bins] = statelevel(y,n)
ymax = max(y);
ymin = min(y)-eps;

% Compute histogram
idx = ceil(n * (y-ymin)/(ymax-ymin));
idx = idx(idx>=1 & idx<=n);
histogram = zeros(n, 1);
for i=1:numel(idx)
    histogram(idx(i)) = histogram(idx(i)) + 1;
end

% Compute center of each bin
ymin = min(y);
Ry = ymax-ymin;
dy = Ry/n;
bins = ymin + ((1:n)-0.5)*dy;

% Compute state levels
iLowerRegion = find(histogram > 0, 1, 'first');
iUpperRegion = find(histogram > 0, 1, 'last');

iLow  = iLowerRegion(1);
iHigh = iUpperRegion(1);

% Define the lower and upper histogram regions halfway between the lowest
% and highest nonzero bins.
lLow  = iLow;
lHigh = iLow + floor((iHigh - iLow)/2);
uLow  = iLow + floor((iHigh - iLow)/2);
uHigh = iHigh;

% Upper and lower histograms
lHist = histogram(lLow:lHigh, 1);
uHist = histogram(uLow:uHigh, 1);

levels = zeros(1,2);
[~, iMax] = max(lHist(2:end));
[~, iMin] = max(uHist);
levels(1) = ymin + dy * (lLow + iMax(1) - 1.5);
levels(2) = ymin + dy * (uLow + iMin(1) - 1.5);
end

function y = nanmean(x,dim)
% FORMAT: Y = NANMEAN(X,DIM)
% 
%    Average or mean value ignoring NaNs
%
%    This function enhances the functionality of NANMEAN as distributed in
%    the MATLAB Statistics Toolbox and is meant as a replacement (hence the
%    identical name).  
%
%    NANMEAN(X,DIM) calculates the mean along any dimension of the N-D
%    array X ignoring NaNs.  If DIM is omitted NANMEAN averages along the
%    first non-singleton dimension of X.
%
%    Similar replacements exist for NANSTD, NANMEDIAN, NANMIN, NANMAX, and
%    NANSUM which are all part of the NaN-suite.
%
%    See also MEAN

% -------------------------------------------------------------------------
%    author:      Jan Glscher
%    affiliation: Neuroimage Nord, University of Hamburg, Germany
%    email:       glaescher@uke.uni-hamburg.de
%    
%    $Revision: 1.1 $ $Date: 2004/07/15 22:42:13 $

if isempty(x)
	y = NaN;
	return
end

if nargin < 2
	dim = min(find(size(x)~=1));
	if isempty(dim)
		dim = 1;
	end
end

% Replace NaNs with zeros.
nans = isnan(x);
x(isnan(x)) = 0; 

% denominator
count = size(x,dim) - sum(nans,dim);

% Protect against a  all NaNs in one dimension
i = find(count==0);
count(i) = ones(size(i));

y = sum(x,dim)./count;
y(i) = i + NaN;
% $Id: nanmean.m,v 1.1 2004/07/15 22:42:13 glaescher Exp glaescher $
end

function y = nanstd(x,dim,flag)
% FORMAT: Y = NANSTD(X,DIM,FLAG)
% 
%    Standard deviation ignoring NaNs
%
%    This function enhances the functionality of NANSTD as distributed in
%    the MATLAB Statistics Toolbox and is meant as a replacement (hence the
%    identical name).  
%
%    NANSTD(X,DIM) calculates the standard deviation along any dimension of
%    the N-D array X ignoring NaNs.  
%
%    NANSTD(X,DIM,0) normalizes by (N-1) where N is SIZE(X,DIM).  This make
%    NANSTD(X,DIM).^2 the best unbiased estimate of the variance if X is
%    a sample of a normal distribution. If omitted FLAG is set to zero.
%    
%    NANSTD(X,DIM,1) normalizes by N and produces the square root of the
%    second moment of the sample about the mean.
%
%    If DIM is omitted NANSTD calculates the standard deviation along first
%    non-singleton dimension of X.
%
%    Similar replacements exist for NANMEAN, NANMEDIAN, NANMIN, NANMAX, and
%    NANSUM which are all part of the NaN-suite.
%
%    See also STD

% -------------------------------------------------------------------------
%    author:      Jan Glscher
%    affiliation: Neuroimage Nord, University of Hamburg, Germany
%    email:       glaescher@uke.uni-hamburg.de
%    
%    $Revision: 1.1 $ $Date: 2004/07/15 22:42:15 $

if isempty(x)
	y = NaN;
	return
end

if nargin < 3
	flag = 0;
end

if nargin < 2
	dim = min(find(size(x)~=1));
	if isempty(dim)
		dim = 1; 
	end	  
end


% Find NaNs in x and nanmean(x)
nans = isnan(x);
avg = nanmean(x,dim);

% create array indicating number of element 
% of x in dimension DIM (needed for subtraction of mean)
tile = ones(1,max(ndims(x),dim));
tile(dim) = size(x,dim);

% remove mean
x = x - repmat(avg,tile);

count = size(x,dim) - sum(nans,dim);

% Replace NaNs with zeros.
x(isnan(x)) = 0; 


% Protect against a  all NaNs in one dimension
i = find(count==0);

if flag == 0
	y = sqrt(sum(x.*x,dim)./max(count-1,1));
else
	y = sqrt(sum(x.*x,dim)./max(count,1));
end
y(i) = i + NaN;

% $Id: nanstd.m,v 1.1 2004/07/15 22:42:15 glaescher Exp glaescher $
end

function df = cent_diff_n(f,h,n)
% df = cent_diff_n(f,h,n)
% Computes an n-point central difference of function f with spacing h.
% Returns a vector df of same size as f.
% Input f must be a vector with evenly spaced points.
% Input n must be 3,5,7, or 9.
% All three inputs are required.
%
% Differences for points near the edges are calculated with lower order.
% For example, if n=5 and length(f)=10, then 3-point central differencing is used
% to calculate values at points 2 and 9, 2-point forward differencing is used for
% point 1, 2-point backward differencing is used for point 10, and 5-point central
% differencing is used for points 3-7.
%
% If f contains less than n points, the order will be downgraded to the
% maximum possible.  Ex: if length(f) = 6, n will be downgraded to 5.
%
% Differencing formulae from: http://www.holoborodko.com/pavel/numerical-methods/numerical-derivative/central-differences/
% Accessed 4/10/12.

if nargin < 3
    error('Not enough inputs.  See help documentation.')
end

if ~isscalar(h)
    error('Input h must be a scalar value.')
end

possible_ns = [3,5,7,9];
if ~ismember(n,possible_ns)
    error('Input n must be 3,5,7, or 9')
end

numPts = length(f);
if numPts < n
    newN = max(possible_ns(possible_ns<=numPts));
    warnstr = [num2str(n) '-point differencing was requested,\n'...
        'but input function only has ' num2str(numPts) ' points.\n'...
        'Switching to ' num2str(newN) '-point differencing.'];
    warning(warnstr,'%s')
    n = newN;
end

df_1 = b_diff(f,h);
df_End = f_diff(f,h);

% Calculate 3-point for all
df_3pt = c_diff(f,h,3);

if n >=5
    df_5pt = c_diff(f,h,n);
    % For the 2nd and next-to-last grid point, use 3-point differencing.
    df_2 = df_3pt(1);
    df_Endm1 = df_3pt(end);
end
if n >=7
    df_7pt = c_diff(f,h,7);
    % For the 3nd and 2nd from last grid point, use 5-point differencing.
    df_3 = df_5pt(1);
    df_Endm2 = df_5pt(end);
end
if n>= 9
    df_9pt = c_diff(f,h,9);
    % For the 4nd and 3rd from last grid point, use 7-point differencing.
    df_4 = df_7pt(1);
    df_Endm3 = df_7pt(end);
end

switch n
    case 3
        df = [df_1 df_3pt df_End];
    case 5
        df = [df_1 df_2 df_5pt df_Endm1 df_End];
    case 7
        df = [df_1 df_2 df_3 df_7pt df_Endm2 df_Endm1 df_End];
    case 9
        df = [df_1 df_2 df_3 df_4 df_9pt df_Endm3 df_Endm2 df_Endm1 df_End];
end
end

function df = c_diff(f,h,n)
midStartPoint = ceil(n/2); % First point at which full n points can be used
midEndPoint = length(f)-midStartPoint+1; % Last point at which full n points can be used

df = [];
for k = midStartPoint:midEndPoint
    switch n
        case 3
            df_k = (f(k+1) - f(k-1))/(2*h);
        case 5
            df_k = (f(k-2) - 8*f(k-1) + 8*f(k+1) - f(k+2))/(12*h);
        case 7
            df_k = (-f(k-3) + 9*f(k-2) - 45*f(k-1) + 45*f(k+1) - 9*f(k+2) + f(k+3))/(60*h);
        case 9
            df_k = (3*f(k-4) - 32*f(k-3) + 168*f(k-2) - 672*f(k-1) + 672*f(k+1) - 168*f(k+2) + 32*f(k+3) - 3*f(k+4))/(840*h);
    end
    df = [df df_k];
end
end

function df1 = b_diff(f,h)
df1 = (f(2)-f(1))/h;
end

function dfEnd = f_diff(f,h)
dfEnd = (f(end)-f(end-1))/h;
end