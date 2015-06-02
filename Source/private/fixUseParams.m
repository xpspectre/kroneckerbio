function [UseParams, nTk] = fixUseParams(UseParams, nk, nCon)
%fixUseParams A helper function for functions that support variently active
%   kinetic parameters to standardize them in a column vector of logicals.
%
%   [UseParams, nTk] = fixUseParams(UseParams, nk)
%
%   Inputs
%   UseParams: [ logical matrix nk by nCon | logical vector nk | positive integer vector ]
%       Indicates the kinetic parameters that will be allowed to vary
%       during the optimization.
%       1) logical matrix nk by nCon to indicate active parameters for each
%           condition
%       2) logical index vector nk to indicate all conditions have the same
%           active parameters
%       3) positive integer vector into nk, all conditions have the same active
%           parameters (TODO: either deprecate or make integer matrix analog)
%   nk: [ nonegative integer scalar ]
%       The number of kinetic parameters in the model
%   nCon: [ nonnegative integer {1} ]
%       Number of experimental conditions
%
%   Outputs
%   UseParams: [ logical matrix nk by nCon ]
%       Standard form of UseParams
%   nTk: [ nonnegative integer scalar ]
%       Number of active kinetic parameters

% (c) 2015 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

if nargin < 3
    nCon = 1;
end

if isnumeric(UseParams)
    if isscalar(UseParams) && isnan(UseParams) % Default
        UseParams = true(nk,nCon);
    else % integer vector
        UseParams = vec(UseParams);
        assert(all(UseParams <= nk), 'KroneckerBio:UseParams:LinearIndexOutOfRange', 'UseParams, when a linear index, can have no value larger than nk. Use a logical matrix to refer to different parameters on each condition.')
        temp = false(nk,nCon);
        temp(UseParams,:) = true;
        UseParams = temp;
    end
elseif islogical(UseParams)
    if length(UseParams) == length(vec(UseParams)) % logical vector
        assert(length(UseParams) == nk, 'KroneckerBio:UseParams:InvalidLogicalVectorSize', 'UseParams, when a logical vector, must have a number of elements equal to nk')
        UseParams = repmat(vec(UseParams), 1, nCon);
    else % logical matrix
        assert(numel(UseParams) == nk*nCon, 'KroneckerBio:UseParams:InvalidLogicalMatrixSize', 'UseParams, when a logical matrix, must have a number of elements equal to nk*nCon')
    end
else
    error('KroneckerBio:UseParams:InvalidType', 'UseParams must be provided as logical or linear index into m.k')
end
nTk = nnz(UseParams);
