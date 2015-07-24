function m = AddOmega(m, names, Omega)
% Add matrix Omega of eta covariances to Model.Analytic. These should match the
% specified etas.
%
% Inputs:
%   m: [ Model.Analytic struct ]
%       The model to which Omega will be added
%   names: [ nEta x 1 cell vector of strings ]
%       Names of parameters that specifies parameters and order of covariance
%       matrix
%   Omega [ nEta x nEta lower triangular double matrix ]
%       Values/initial guesses of variance of etas. Specify as lower triangular
%       matrix. If not, warn and only take lower triangular elements.
% Outputs:
%   m: [ Model.Analytic struct ]
%       The model with Omega added

% Clean up inputs
names = vec(names);
n = length(names);
if nargin < 3
    Omega = zeros(n);
end

% TODO: Sanity checks:
%   Verify that names are present in model and have corresponding etas
%   Vefity that Omega is lower triangular. If not, warn and only take lower triangular elements.

% Add sigmas going down cols of lower triangular Omega
for j = 1:n
    for i = j:n
        nk = m.add.nk + 1;
        m.add.nk = nk;
        m.add.Parameters = growParametersAnalytic(m.add.Parameters, nk);
        
        m.add.Parameters(nk).Name = fixParameterName(['omega__' names{i} '__' names{j}]);
        m.add.Parameters(nk).ID = genUID;
        m.add.Parameters(nk).Value = Omega(i,j);
%         fprintf([num2str(i) ',' num2str(j) '\n'])
    end
end

m.Ready = false;