function [esttype, arglist] = psdesttype(validtypes,arglist)
%PSDESTTYPE - return the PSD estimation type option and remove it from the
%argument list
%
%   validtypes  - a cell array of valid estimator types
%                 (e.g. {'power','ms','psd'})
%   varargin    - the input argument list
%
%   Errors out if different estimation types are specified in varargin.    

%   Copyright 2012 The MathWorks, Inc.

% use 'psd' by default
esttype = 'psd'; 
found = false; 

for i=1:numel(validtypes) 
matches = find(strcmpi(validtypes{i},arglist)); 
if ~isempty(matches) 
if ~found
    found = true;
    esttype = validtypes{i};
    arglist(matches) = [];
    else
    error(message('signal:psdoptions:ConflictingEstTypes', ...
                  esttype,validtypes{i}));
    end
end
end 