function m = AddSigma(m, names, Sigma)
% Add matrix Sigma of eps covariances to Model.Analytic. These should match the
% specified eps.
%
% Inputs:
%   m: [ Model.Analytic struct ]
%       The model to which Sigma will be added
%   names: [ nEps x 1 cell vector of strings ]
%       Names of parameters that specifies parameters and order of covariance
%       matrix
%   Sigma [ nEps x nEps lower triangular double matrix ]
%       Values/initial guesses of variance of etas. Specify as lower triangular
%       matrix. If not, warn and only take lower triangular elements.
% Outputs:
%   m: [ Model.Analytic struct ]
%       The model with Sigma added

% Clean up inputs
names = vec(names);
n = length(names);
if nargin < 3
    Sigma = zeros(n);
end

% TODO: Sanity checks:
%   Verify that names are present in model and have corresponding etas
%   Vefity that Sigma is lower triangular. If not, warn and only take lower triangular elements.

% Add sigmas going down cols of lower triangular Sigma
for j = 1:n
    for i = j:n
        m = AddParameter(m, ['sigma__' names{i} '__' names{j}], Sigma(i,j));
    end
end

m.Ready = false;