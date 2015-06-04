function [UseParams, nTk] = fixUseParams(UseParams, nk, nCon)
%fixUseParams A helper function for functions that support variently active
%   kinetic parameters to standardize them in a column vector of logicals.
%
%   [UseParams, nTk] = fixUseParams(UseParams, nk, nCon)
%
%   Inputs
%   UseParams [ integer matrix nk by nCon | logical matrix nk by nCon | logical vector nk ]
%       Indicates the kinetic parameters that will be allowed to vary
%       during the optimization. UseParams can be:
%       1) integer matrix nk by nCon to indicate active parameters for each
%           condition. In each row, different integers indicate rate
%           parameters that differ between conditions. 0 entries indicate
%           not fit. [Standard form]
%       2) logical matrix nk by nCon to indicate active parameters for each
%           condition, with all parameters unique.
%       3) logical index vector nk to indicate all conditions (or single condition) have the same
%           active parameters, with all parameters unique
%   nk [ nonegative integer scalar ]
%       The number of kinetic parameters in the model
%   nCon [ nonnegative integer {1} ]
%       Number of experimental conditions
%
%   Outputs
%   UseParams [ integer matrix nk by nCon | integer vector nk ]
%       Standard form of UseParams. Matrix indicates parameters allowed to vary
%       between conditions according to unique integers. Vector indicates single
%       parameter set.
%   nTk [ nonnegative integer scalar ]
%       Number of active kinetic parameters

% (c) 2015 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

if isnumeric(UseParams)
    if isscalar(UseParams) && isnan(UseParams) % Default
        UseParams = true(nk,1);
    else % integer matrix
        assert(all(size(UseParams) == [nk,nCon]), 'KroneckerBio:UseParams:InvalidIntegerMatrixSize', 'UseParams, when an integer matrix, must have dimension nk by nCon')
    end
elseif islogical(UseParams)
    if length(UseParams) == length(vec(UseParams)) % logical vector
        assert(length(UseParams) == nk, 'KroneckerBio:UseParams:InvalidLogicalVectorSize', 'UseParams, when a logical vector, must have a number of elements equal to nk')
        temp = vec(1:nk);
        temp(~UseParams) = 0;
        UseParams = temp;
    else % logical matrix
        assert(all(size(UseParams) == [nk,nCon]), 'KroneckerBio:UseParams:InvalidLogicalMatrixSize', 'UseParams, when a logical matrix, must have dimensions nk by nCon')
        temp = reshape(1:nk*nCon, nk, nCon);
        temp(~UseParams) = 0;
        UseParams = temp;
    end
else
    error('KroneckerBio:UseParams:InvalidType', 'UseParams must be provided as logical or linear index into m.k')
end
nTk = getnTk(UseParams);
