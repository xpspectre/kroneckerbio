function [UseSeeds, nTs] = fixUseSeeds(UseSeeds, ns, nCon)
%fixUseSeeds A helper function for functions that support variant initial
%   concentrations in experiments. It converts a variety of ways to
%   specific the active seed parameters and standardizes them into
%   matrix of logicals indicating the active parameters.
%
%   [UseSeeds, nTs] = fixUseSeeds(UseSeeds, UseModelSeeds, ns, nCon)
%
%   Inputs
%   UseSeeds: [ logical matrix ns by nCon | logical vector ns | positive integer vector ]
%       Indicates the seed parameters that will be allowed to vary
%       during the optimization.
%       1) logical matrix ns by nCon to indicate active seed parameters for each
%           condition
%       2) logical index vector ns to indicate all conditions have the same
%           active seed parameters
%       3) positive integer vector into ns, all conditions have the same active
%           seed parameters (TODO: either deprecate or make integer matrix analog)
%   ns: [ nonegative integer scalar ]
%       The number of seed parameters in the model
%   nCon: [ nonnegative integer {1} ]
%       Number of experimental conditions
%
%   Outputs
%   UseSeeds: [ logical matrix ns by nCon ]
%       Standard form of UseSeeds
%   nTs: [ nonnegative integer ]
%       Number of active seed parameters

% (c) 2015 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

if nargin < 3
    nCon = 1;
end

if isnumeric(UseSeeds)
    if isscalar(UseSeeds) && isnan(UseSeeds) % Default
        UseSeeds = true(ns,nCon);
    else % integer vector
        UseSeeds = vec(UseSeeds);
        assert(all(UseSeeds <= ns), 'KroneckerBio:UseSeeds:LinearIndexOutOfRange', 'UseSeeds, when a linear index, can have no value larger than ns. Use a logical matrix to refer to different parameters on each condition.')
        temp = false(ns,nCon);
        temp(UseSeeds,:) = true;
        UseSeeds = temp;
    end
elseif islogical(UseSeeds)
    if length(UseSeeds) == length(vec(UseSeeds)) % logical vector
        assert(length(UseSeeds) == ns, 'KroneckerBio:UseSeeds:InvalidLogicalVectorSize', 'UseSeeds, when a logical vector, must have a number of elements equal to ns')
        UseSeeds = repmat(vec(UseSeeds), 1, nCon);
    else % logical matrix
        assert(numel(UseSeeds) == ns*nCon, 'KroneckerBio:UseSeeds:InvalidLogicalMatrixSize', 'UseSeeds, when a logical matrix, must have a number of elements equal to ns*nCon')
    end
else
    error('KroneckerBio:UseSeeds:InvalidType', 'UseSeeds must be provided as logical or linear index into con.s')
end
nTs = nnz(UseSeeds);
