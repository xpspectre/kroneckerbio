function obs = observationFocei(outputs, times, method, name)
% Nonlinear mixed effects observation and objective function for an
% individual's theta
%
% Inputs:
%   outputs [ string | nOutputs cell array of strings ]
%       Name of output{s} of interest
%   times [ nt x 1 vector of nonnegative doubles ]
%       Vector of times of interest
%   method [ {'FO'} | 'FOCEI' ]
%       Approximation method
%   name [ string {'Obs_Individual'} ]
%       String name of model
%
% Outputs:
%   obj [ Observation.Data.Focei struct ]
%       Observation struct
%
% Notes:
%   Although it supports multiple outputs/measurements per individual, they have
%   to be at the same times (for now)

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
    name = 'Obs_Individual';
end

% Validate method
switch method
    case 'FO'
        % pass
    case 'FOCEI'
        % pass
    otherwise
        error('KroneckerBio:observationFocei:InvalidMethod', 'Method %s is invalid. Use FO or FOCEI.', method)
end

% Clean up outputs
outputs = fixOutputs(outputs);

% Find unique timelist
times = row(unique(times));

% Sanity checks
nOutputs = length(outputs);
nt = length(times);
assert(iscellstr(outputs), 'KroneckerBio:observationFocei:OutputsNotCellString', 'Outputs is not a cell vector of strings')
assert(all(size(outputs) == [nOutputs,1]), 'KroneckerBio:observationFocei:InvalidOutputsFormat', 'Outputs is not a column vector')
assert(all(size(times) == [1,nt]), 'KroneckerBio:observationFocei:InvalidTimesFormat', 'Times is not a row vector')

obs = [];
obs.Type    = 'Observation.Data.Focei';
obs.Name    = name;
obs.Complex = false;

obs.tF            = max([0, times]);
obs.DiscreteTimes = times;
obs.Outputs       = outputs;

obs.Simulation = @simulation;
obs.Objective  = @objective;

obs = pastestruct(observationZero(), obs);

    function sim = simulation(int)
        sim.Type   = 'Simulation.Data.Focei';
        sim.Name   = name;
        sim.Output = outputs;
        sim.Times  = times;
        sim.Solution = int;
    end

    function obj = objective(measurements)
        obj = objectiveFocei(outputs, times, method, name, measurements);
    end

end

function obj = objectiveFocei(outputs, times, method, name, measurements)
% Objective function corresponding to observation
%
% Inputs:
%   outputs [ nOutputs x 1 cell vector of strings ]
%   times [ 1 x ni double vector ]
%   measurements [ nOutputs x ni double matrix ]
%   method [ 'FO' | 'FOCEI' ]
%   name [ string ]
%
% Outputs:
%   obj [ Objective.Data.Focei struct ]

% Verify measurements match times
nOutputs = length(outputs);
nt       = length(times);
assert(all(size(measurements) == [nOutputs, nt]), 'KroneckerBio:objectiveFocei:InvalidMeasurements', 'Dimensions of measurements don''t match nTimes and/or nOutputs')

% Inherit observation
obj = observationFocei(outputs, times, method, name);

obj.Type = 'Objective.Data.Focei';
obj.Continuous = false;

obj.G      = @G;
obj.dGdk   = @dGdk;
obj.d2Gdk2 = @d2Gdk2;

obj.p    = @p;
obj.logp = @logp;
obj.F    = @F;
obj.Fn   = @Fn;

private = [];
private.outputs      = outputs;
private.times        = times;
private.measurements = measurements;
obj.private = private;

obj = pastestruct(objectiveZero(), obj);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outer parameter fitting functions %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [val, discrete] = G(int)
        % Calculate outer objective function value
        %         fprintf('Calculating obj fun with 1st order sensitivities\n')
        order = 1;
        layer = 'outer';
        val = calcFocei(int, outputs, measurements, order, layer);
        discrete = 0;
    end

    function val = dGdk(int)
        % Calculate outer objective function gradient w.r.t. theta
        % Note: Since NLME requires that the n+1 order sensitivities are calculated
        %   for each n order calculation (i.e., 1st order sensitivities for G
        %   calculation), this will be called extraneously in the current scheme
        %   when only G is needed. Detect that the calling function wants only a
        %   the lower order term by looking for order of sensitivities.
        % TODO: fix the above point so this is called exactly when it needs.
        %         fprintf('Calculating obj fun and grad with 2nd order sensitivities\n')
        order = 2;
        layer = 'outer';
        thetaInds = getParameterInds(int.k_names);
        val = zeros(int.nk, 1);
        if isfield(int, 'd2ydx2')
            [~, dGdtheta] = calcFocei(int, outputs, measurements, order, layer);
            val(thetaInds) = dGdtheta;
        end
    end

    function val = d2Gdk2(int)
        % Calculate outer objective function Hessian w.r.t. theta
        % Note: This is not needed for fitting and hasn't been derived analytically.
        %   It is called extraneously when dGdk is desired, so return a dummy
        %   matrix of 0's to not add anything.
        order = 3;
        layer = 'outer';
        val = zeros(int.nk);
        if isfield(int, 'd3ydx3') % calling function wants Hessian - will never be called
            [~, ~, d2Gdtheta2] = calcFocei(int, outputs, measurements, order, layer);
        end
    end

end
