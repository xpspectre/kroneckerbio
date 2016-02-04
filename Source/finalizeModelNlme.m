function m = finalizeModelNlme(m, varargin)
%finalizeModelNlme prepares a constructed Model.Nlme for use with the
% rest of kroneckerbio and nonlinear mixed effecst modeling.
% Generates symbolic expressions and requires the symbolic
% toolbox.

% Sanity checks
assert(isfield(m, 'OmegaNames'), 'kroneckerbio:finalizeModelNlme:MissingOmegaNamesField')

% Specific Model.Nlme finalization
% Symbolically calculate Omega, inverses, and derivative and turn into
%   function handles
kNames = vec({m.Parameters.Name});
fOmega = calcOmega(m.OmegaNames, kNames);

% General Model.Analytic finalization
m = finalizeModelAnalytic(m, varargin{:});

% Augment finalized model with Nlme-specific fields
m.fOmega = fOmega;