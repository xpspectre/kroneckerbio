function obj = objectiveNLME(output, times, measurements, name)
% Dummy objective function for nonlinear mixed effects fitting. Not composable
% with regular objective functions because this isn't meant to be used with
% FitObjective.
%
% Inputs:
%   output [ string ]
%       Name of output of interest
%   times [ nt x 1 vector of nonnegative doubles ]
%       Vector of times of interest
%   measurements [ nt x 1 vector of nonnegative doubles ]
%       Vector of measurements of interest
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
end

if isempty(name)
    name = 'ObjectiveNLME';
end

% Make sure a single output name is provided
assert(ischar(output), 'objectiveNLME:invalidOutput', 'Output must be a string name of the output')

% Make sure times line up with measurements
times = vec(times);
measurements = vec(measurements);
assert(length(times) == length(measurements), 'objectiveNLME:measurementsMismatch', 'Number of times and measurements don''t match')

obj.Type = 'Dummy.Objective.NLME';
obj.Name = name;

obj.Outputs = output;
obj.Times = times;
obj.Measurements = measurements;

end

