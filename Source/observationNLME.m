function obs = observationNLME(output, timelist, method, name)
% Nonlinear mixed effects observation and objective function for a patient
%
% Inputs:
%   output [ string ]
%       Name of output of interest
%   timelist [ nt x 1 vector of nonnegative doubles ]
%       Vector of times of interest
%   method [ {'FO'} | 'FOCE' ]
%       Approximation method
%   name [ string {'NLME'} ]
%       String name of model
% Outputs:
%   obj [ Objective.NLME struct ]
%       Observation struct
%
% Note/TODO: Currently only supports a single output. Extend to use an arbitrary
% number of outputs later.

% Clean up inputs
if nargin < 4
    name = [];
    if nargin < 3
        method = [];
    end
end

if isempty(method)
    method = 'FO';
end
if isempty(name)
    name = 'NLME';
end

% Make sure a single output name is provided
assert(ischar(output), 'observationNLME_FO:nonscalarOutput', 'Output must be a string')

% Find unique timelist
discrete_times = row(unique(timelist));

obs = [];
obs.Type    = 'Observation.Data.NLME';
obs.Name    = name;
obs.Complex = false;

obs.tF            = max([0, discrete_times]);
obs.DiscreteTimes = discrete_times;
obs.Outputs       = output;

obs.Simulation = @simulation;
obs.Objective  = @objective;

obs = pastestruct(observationZero(), obs);

    function sim = simulation(int)
        sim.Type = 'Simulation.Data.NLME';
        sim.Name = name;
        sim.Output = output;
        sim.Times = discrete_times;
        
        y_Ind = lookupOutput(output, int); % take actual output being compared to measurements (w/ error model from eps)
        sim.Predictions = int.y(y_Ind,:);
    end

    function obj = objective(measurements)
        obj = objectiveNLME(output, timelist, measurements, method, name);
    end

end

function obj = objectiveNLME(output, timelist, measurements, method, name)

% Make sure a single output name is provided
assert(ischar(output), 'observationNLME_FO:nonscalarOutput', 'Output must be a string')

% Find unique timelist
discrete_times = row(unique(timelist));

% Inherit observation
obj = observationNLME(output, timelist, method, name);

obj.Type = 'Objective.Data.NLME';
obj.Continuous = false;

obj.G      = @G;
obj.dGdk   = @dGdk;
obj.d2Gdk2 = @d2Gdk2;

obj.p    = @p;
obj.logp = @logp;
obj.F    = @F;
obj.Fn   = @Fn;

obj = pastestruct(objectiveZero(), obj);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Parameter fitting functions %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function [val, discrete] = G(int)
        
        % Get output Fi__{output}
        [y_Ind, Fi_y_Ind] = lookupOutput(output, int);
        ny = int.ny;
        Fi = int.y(Fi_y_Ind,:)';
        
        % Get objective function value for each patient
        RES = measurements - Fi;
        
        % Extract indexing for NLME fit
        kNames = int.k_names;
        k = int.k;
        [thetaInds, etaInds, omegaInds, epsInds, sigmaInds] = getNLMEInds(kNames);
        
        OM = getOmega(k, kNames, thetaInds, omegaInds);
        SG = getSigma(k, kNames, epsInds, sigmaInds);
        
        Gi = getGi(int.dydT, kNames, thetaInds, etaInds, Fi_y_Ind, ny);
        Hi = getHi(int.dydT, epsInds, y_Ind, ny);
        
        COV = Gi*OM*Gi' + diag(diag(Hi*SG*Hi'));
        WRES = matSqrtInv(COV)*RES;
        
        val = log(det(COV)) + WRES'*WRES;
        
        % Return discrete times as well
        discrete = discrete_times;
    end

    function val = dGdk(int)
        warning('Analytic derivative for NLME dGdk not implemented yet. Derive and implement it or approximate with finite differences.')
        val = zeros(int.nk, 1);
    end

    function val = d2Gdk2(int)
        warning('Analytic derivative for NLME d2Gdk2 not implemented yet. Derive and implement it or approximate with finite differences.')
        val = zeros(int.nk, int.nk);
    end


end

%% Helper functions
function [y_Ind, Fi_y_Ind] = lookupOutput(name, int)
y_Ind    = find(ismember(int.y_names, name));
Fi_y_Ind = find(ismember(int.y_names, ['Fi__' name]));
end

function [thetaInds, etaInds, omegaInds, epsInds, sigmaInds] = getNLMEInds(kNames)
% Get indices of NLME components from rate params list
nk = length(kNames);
thetaInds = zeros(nk,1);
etaInds   = zeros(nk,1);
omegaInds = zeros(nk,1);
epsInds   = zeros(nk,1);
sigmaInds = zeros(nk,1);
for i = 1:nk
    kName = kNames{i};
    if ~isempty(regexp(kName, '^eta__', 'once'))
        etaInds(i) = i;
    elseif ~isempty(regexp(kName, '^omega__', 'once'))
       omegaInds(i) = i;
    elseif ~isempty(regexp(kName, '^eps__', 'once'))
       epsInds(i) = i;
    elseif ~isempty(regexp(kName, '^sigma__', 'once'))
       sigmaInds(i) = i;
    else % not a NLME param, so call a theta
        thetaInds(i) = i;
    end
end
thetaInds(thetaInds==0) = [];
etaInds(etaInds==0)     = [];
omegaInds(omegaInds==0) = [];
epsInds(epsInds==0)     = [];
sigmaInds(sigmaInds==0) = [];
end

function Omega = getOmega(k, kNames, thetaInds, omegaInds)
% Assemble full Omega matrix by looking up params from k and kNames. Omega is a
% symmetric n x n matrix where n = number of "eta" parameters being fit and
% order the same as k.
% Parameters w/o inter-individual variability eta specified have omegas set = 0.

omegas = k(omegaInds);
omegaNames = kNames(omegaInds);
thetaNames = kNames(thetaInds);

no = length(omegaInds);
n = length(thetaInds);
Omega = zeros(n);
for i = 1:no
    tokens = regexp(omegaNames{i}, '^omega__(.+)__(.+)$', 'tokens', 'once')';
    assert(iscell(tokens) && length(tokens) == 2, 'observationNLME:getOmega:invalidOmegaTokens', 'Parser didn''t find correct params omega name %s is referring to', omegaNames{i})
    ix = find(ismember(thetaNames, tokens{1}));
    jx = find(ismember(thetaNames, tokens{2}));
    
    Omega(ix,jx) = omegas(i);
    Omega(jx,ix) = omegas(i);
end
end

function Sigma = getSigma(k, kNames, epsInds, sigmaInds)
% Assemble full Sigma matrix by looking up params from k and kNames. Sigma is a
% symmetric n x n matrix where n = number of "eps" parameters being fit.

sigmas = k(sigmaInds);
sigmaNames = kNames(sigmaInds);
epsNames   = kNames(epsInds);

% Trim leading "eps__" from epsNames
n = length(epsInds);
for i = 1:n
    epsNames{i} = strrep(epsNames{i}, 'eps__', '');
end

no = length(sigmaInds);
Sigma = zeros(n);
for i = 1:no
    tokens = regexp(sigmaNames{i}, '^sigma__(.+)__(.+)$', 'tokens', 'once')';
    assert(iscell(tokens) && length(tokens) == 2, 'observationNLME:getSigma:invalidSigmaTokens', 'Parser didn''t find correct params sigma name %s is referring to', sigmaNames{i})
    ix = find(ismember(epsNames, tokens{1}));
    jx = find(ismember(epsNames, tokens{2}));
    
    Sigma(ix,jx) = sigmas(i);
    Sigma(jx,ix) = sigmas(i);
end
end

function Gi = getGi(dydT, kNames, thetaInds, etaInds, Fi_y_Ind, ny)
% Gi = dF/deta, a nTimes (nObservations) x netas matrix
% dydT represents the ny x nT x nt matrix of partial derivatives of
%   outputs wrt all parameters being optimized over as (ny x nT) * nt. Going
%   down dydT, y changes fastest, i.e., the first ny entries are dy1/dT1, dy2/dT1, dy3/dT1, etc.
%   the next ny entries are dy1/dT2, etc.
% kNames is sufficient because k's are first in T and other types of params aren't used

thetaNames = kNames(thetaInds);
etaNames   = kNames(etaInds);

nTheta = length(thetaInds);
nEta = length(etaInds);

% Get ks corresponding to etas
dydT_eta_Ind = zeros(nTheta,1);
for i = 1:nEta
    token = regexp(etaNames{i}, '^eta__(.+)$', 'tokens', 'once');
    thetaNameThisEta = token{1};
    thetaMaskThisEta = ismember(thetaNames, thetaNameThisEta);
    dydT_eta_Ind(thetaMaskThisEta) = find(ismember(kNames, etaNames{i}));
end

% Assemble Gi matrix
% Thetas with no corresponding etas will have 0 values
nt = size(dydT,2);
Gi = zeros(nt,nTheta);
for i = 1:nTheta
    offset = ny*(dydT_eta_Ind(i)-1); % to find the right parameter block
    ind = Fi_y_Ind + offset; % to find the right output
    if ind ~= 0
        Gi(:,i) = dydT(ind,:)';
    end
end
end

function Hi = getHi(dydT, epsInds, y_Ind, ny)
% Hi = dY/deps
nEps = length(epsInds);
nt = size(dydT,2);
Hi = zeros(nt,nEps);
for i = 1:nEps
    offset = ny*(epsInds(i)-1); % to find the right parameter block
    ind = y_Ind + offset; % to find the right output
    Hi(:,i) = dydT(ind,:)';
end
end
