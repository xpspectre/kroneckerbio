function [UseParams, nTk] = fixUseParams(UseParams, nk)
%fixUseParams A helper function for functions that support variently active
%   kinetic parameters to standardize them in a column vector of logicals.
%
%   [UseParams, nTk] = fixUseParams(UseParams, nk)
%
%   Inputs
%   UseParams: [ logical vector nk | positive integer vector ]
%       Indicates the kinetic parameters that will be allowed to vary
%   nk: [ nonegative integer scalar ]
%       The number of kinetic parameters in the model
%
%   Outputs
%   UseParams: [ logical vector nk ]
%       Standard form of UseParams
%   nTk: [ nonnegative integer scalar ]
%       Number of active kinetic parameters

% (c) 2015 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

if islogical(UseParams)
    UseParams = vec(UseParams);
elseif isnumeric(UseParams)
    temp = UseParams;
    UseParams = zeros(nk, 1);
    UseParams(temp) = 1;
    UseParams = logical(UseParams);
else
    error('KroneckerBio:fixUseParams', 'UseParams must be a vector of logicals or numerics')
end
nTk = sum(UseParams);
