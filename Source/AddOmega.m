function m = AddOmega(m, names, Omega)
% Add matrix Omega of Eta covariances to Model.Analytic. These should match the
% specified etas. Must be run after thetas and etas are added. Omega is made of
% omegas^2, the variances/covariances.
%
% Inputs:
%   m: [ Model.Analytic struct ]
%       The model to which Omega will be added
%   names: [ nEta x 1 cell vector of strings ]
%       Names of parameters that specifies parameters and order of covariance
%       matrix
%   Omega [ nEta x nEta lower triangular double matrix {10% of corresponding theta values} ]
%       Values/initial guesses of variance of etas. Given as std devs omegas.
%       Specify as lower triangular matrix. If not, warn and only take lower
%       triangular elements. If not specified, fills in default values equaling
%       10% of the corresponding parameter theta values from the model.
%
% Outputs:
%   m: [ Model.Analytic struct ]
%       The model with Omega added

% Clean up inputs
names = vec(names);
n = length(names);
if nargin < 3
    Omega = [];
end

% Verify that parameters are present in model and have corresponding etas
k_names = vec({m.Parameters.Name});
k_names = k_names(~cellfun('isempty', k_names));
[foundNameMask, foundNameInds] = ismember(names, k_names);
assert(all(foundNameMask), 'AddOmega:InvalidName', 'Names: %s not found in model', cellstr2str(names(~foundNameMask)))

% Set default values of variances if not given
if isempty(Omega)
    warning('AddOmega:OmegaNotSpecified', 'Omega values not specified. Setting initial guess of variances to 10%% of Theta values and covariance to 0.')
    k_vals = vec([m.Parameters.Value]); % automatically discards empty entries
    foundVals = k_vals(foundNameInds);
    Omega = diag((0.1*foundVals).^2);
end

% Verify that number of names = dim of Omega
assert(size(Omega,1) == size(Omega,2) && size(Omega,1) == numel(names), 'AddOmega:NumNamesOmegaMismatch', 'Number of names = %g doesn''t match dims of Omega = %g', numel(names), size(Omega,1))

% Verify that Omega is lower triangular. If not, warn and only take lower triangular elements.
if ~istril(Omega)
    warning('AddOmega:NotLowerTriangular', 'Specified Omega was not lower triangular. Only taking lower triangular elements')
    Omega = tril(Omega);
end

% Verify that Omega is a valid covariance matrix (positive semidefinite)
%   TODO: make sure this always holds during the optimization
OmegaFull = tril(Omega,-1)' + Omega; % expand to full symmetric matrix
if any(eig(OmegaFull) < 0)
    warning('AddOmega:OmegaNotPositiveSemiDefinite', 'The full Omega is not positive semi-definite and is an invalid covariance matrix. Using anyway')
end

% Add omegas as individual parameters
for i = 1:n
    for j = 1:i
        m = AddParameter(m, ['omega__' names{i} '__' names{j}], sqrt(Omega(i,j)));
    end
end

m.Ready = false;