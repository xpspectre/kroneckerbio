function obs = observationIndividualVariability(outputs, times, name)
% Nonlinear mixed effects observation and objective function for an
% individual's eta
%
% Inputs:
%   outputs [ string | nOutputs cell array of strings ]
%       Name of output{s} of interest
%   times [ nt x 1 vector of nonnegative doubles ]
%       Vector of times of interest
%   name [ string {'Obs_Individual_Variability'} ]
%       String name of model
%
% Outputs:
%   obj [ Observation.Data.Focei.IndividualVariability struct ]
%       Observation struct
%
% Notes:
%   Although it supports multiple outputs/measurements per individual, they have
%   to be at the same times (for now)

% Clean up inputs
if nargin < 4
    name = [];
end
if isempty(name)
    name = 'Obs_Individual_Variability';
end

% Clean up outputs
outputs = fixOutputs(outputs);

% Find unique timelist
times = row(unique(times));

% Sanity checks
nOutputs = length(outputs);
nt = length(times);
assert(iscellstr(outputs), 'KroneckerBio:observationIndividualVariability:OutputsNotCellString', 'Outputs is not a cell vector of strings')
assert(all(size(outputs) == [nOutputs,1]), 'KroneckerBio:observationIndividualVariability:InvalidOutputsFormat', 'Outputs is not a column vector')
assert(all(size(times) == [1,nt]), 'KroneckerBio:observationIndividualVariability:InvalidTimesFormat', 'Times is not a row vector')

obs = [];
obs.Type    = 'Observation.Data.Focei.IndividualVariability';
obs.Name    = name;
obs.Complex = false;

obs.tF            = max([0, times]);
obs.DiscreteTimes = times;
obs.Outputs       = outputs;

obs.Simulation = @simulation;
obs.Objective  = @objective;

obs = pastestruct(observationZero(), obs);

    function sim = simulation(int)
        sim.Type   = 'Simulation.Data.Focei.IndividualVariability';
        sim.Name   = name;
        sim.Output = outputs;
        sim.Times  = times;
        sim.Solution = int;
    end

    function obj = objective(measurements)
        obj = objectiveIndividualVariability(outputs, times, name, measurements);
    end

end

function obj = objectiveIndividualVariability(outputs, times, name, measurements)
% Objective function for eta optimization for each individual
%
% Inputs:
%   outputs [ nOutputs x 1 cell vector of strings ]
%   times [ 1 x ni double vector ]
%   measurements [ nOutputs x ni double matrix ]
%   name [ string ]
%
% Outputs:
%   obj [ Objective.Data.Focei.IndividualVariability struct ]

% Verify measurements match times
nOutputs = length(outputs);
nt       = length(times);
assert(all(size(measurements) == [nOutputs, nt]), 'KroneckerBio:objectiveIndividualVariability:InvalidMeasurements', 'Dimensions of measurements don''t match nTimes and/or nOutputs')

% Inherit observation
obj = observationIndividualVariability(outputs, times, name);

obj.Type = 'Objective.Data.Focei.IndividualVariability';
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
% Inner parameter fitting functions %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [val, discrete] = G(int)
        % Calculate inner objective function value of optimizing eta
        %         fprintf('Calculating obj fun with 1st order sensitivities\n')
        order = 1;
        layer = 'inner';
        val = calcFocei(int, outputs, measurements, order, layer);
        discrete = 0;
    end

    function val = dGdk(int)
        % Calculate inner objective function gradient w.r.t. eta
        %         fprintf('Calculating obj fun and grad with 2nd order sensitivities\n')
        order = 2;
        layer = 'inner';
        [~, etaInds] = getParameterInds(int.k_names);
        val = zeros(int.nk, 1);
        if isfield(int, 'd2ydx2')
            [~, dGdeta] = calcFocei(int, outputs, measurements, order, layer);
            val(etaInds) = dGdeta;
        end
    end

    function val = d2Gdk2(int)
        % Calculate inner objective function Hessian w.r.t. eta
        order = 3;
        layer = 'inner';
        val = zeros(int.nk);
        if isfield(int, 'd3ydx3') % calling function wants Hessian - will never be called
            [~, ~, d2Gdeta2] = calcFocei(int, outputs, measurements, order, layer);
        end
    end

end
