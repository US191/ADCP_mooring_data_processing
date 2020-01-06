function w = hamm(n)
%HAMM HAMM(N) returns the N-point Hamming window.

%   Copyright (c) 1988-96 by The MathWorks, Inc.
%       $Revision: 1.7 $  $Date: 1996/07/25 16:36:49 $

if n > 1
    w = .54 - .46*cos(2*pi*(0:n-1)'/(n-1));
elseif n == 1
    w = .08;
else
    error('N must be greater than 0.')
end

