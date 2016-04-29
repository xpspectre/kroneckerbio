function [mOut, conOut, obj, opts] = randomExperimentFittingData(m, con, opts, nCon, tF, nMeasurements)
% Generate random experiments for testing objective-based Kronecker
%   functions with multiple experiments
% 
% Takes a model-experiment combination (m and con), chooses random k, s, q,
%   and h from those indicated in opts.Use*, and varies them randomly up and
%   down up to twofold to generate one variant model and a set of
%   nExperiments variant experiments. The function then simulates the variant
%   model and experiments up to time tF, gathers data from nMeasurements random outputs and
%   time points, and creates an objective function array obj fitting to the
%   data with added error. Returns the model(s), experimental conditions,
%   objective functions, and options used to generate the data. To run a fit,
%   use m, nCon x con, obj, and opts; and compare to mOut and conOut.
%
% Inputs:
%   m [ scalar model struct ]
%       Base model to generate experiments from
%   con [ scalar condition struct ]
%       Base condition to generate experiments from
%   opts [ scalar options struct ]
%       Options struct containing parameter fitting spec in the "old" format. Fields:
%       .UseParams [ nk x 1 logical vector ]
%           Parameters to fit. All conditions share these.
%       .UseSeeds [ ns x 1 logical vector | ns x nCon logical matrix]
%           Seeds to fit by condition. All seeds specified here are independent
%           values (allowing all conditions to have different seeds).
%       .UseInputControls [ nCon x 1 cell vector of nqi x 1 logical vectors ]
%           Input control parameters to fit, specified for the inputs for each
%           experimental condition.
%       .UseDoseControls [ nCon x 1 cell vector of nhi x 1 logical vectors ]
%           Dose control parameters to fit, specified for the inputs for each
%           experimental condition.
%   nCon [ scalar positive integer ]
%       Number of experimental conditions to generate off of m
%   tF [ scalar positive double ]
%       Final time for simulated data
%   nMeasurements [ scalar positive integer ]
%       Number of timepoints at which to take measurements of simulated data to
%       fit
%
% Outputs:
%   mOut [ scalar model struct ]
%       Model with randomly modified parameters for comparison
%   conOut [ nCon x 1 model struct vector ]
%       Conditions with randomly modified parameters for comparison
%   obj [ nCon x 1 objective function struct vector ]
%       Objective functions for fitting generated data. Makes a single
%       observationLinearWeightedSumOfSquares obj fun for each con.
%   opts [ scalar options struct ]
%
% Warning: AbsTol is set to default value of 1e-9 no matter what to avoid
%   dealing with standardization, which is complicated.
%
% (c) 2016 David Flowers
% This work is released under the MIT license.

% Input checks
assert(isscalar(m), 'm must be a scalar.')
assert(isscalar(con), 'con must be a scalar.')

usefields = {'UseParams','UseSeeds','UseInputControls','UseDoseControls'};
hasfields = isfield(opts, usefields);
assert(all(hasfields), 'getRandomModelFittingData:MissingUseFields',...
    'Fields %s are required in opts but were not found',...
    strjoin(usefields(~hasfields), ', '))

% Get standardization functions from private directory (bad practice here)
try
    currentdir = pwd;
    cd(fullfile('..','Source','private'))
    fixUseParams_ = @fixUseParams;
    fixUseSeeds_ = @fixUseSeeds;
    fixUseControls_ = @fixUseControls;
    cd(currentdir)
catch ME % Avoids getting stuck in a different directory if an error occurs
    cd(currentdir)
    rethrow(ME)
end

% Standardize UseParams as nk x 1 logical vector, UseSeeds as ns x nCon logical matrix for UseSeeds or nCon x 1
% cell array of n(q or h) x 1 logical vectors otherwise
opts.UseParams = fixUseParams_(opts.UseParams, m.nk);
opts.UseSeeds = fixUseSeeds_(opts.UseSeeds, m.ns, numel(con));
opts.UseInputControls = fixUseControls_(opts.UseInputControls, numel(con), [con.nq]);
opts.UseDoseControls = fixUseControls_(opts.UseDoseControls, numel(con), [con.nh]);

% Set scalar AbsTol to avoid the complexity (obviously this isn't a great
% solution)
opts.AbsTol = 1e-9;

% Function that selects Use* parameters at random and generates a new Use*
% logical vector
selectUses = @(use) use & (randi([0,1],size(use,1),1) == 1);

% Assign random model parameters to be fit
opts.UseParams = selectUses(opts.UseParams);

% Choose parameters to randomize and fit
UseNames = {'UseSeeds','UseInputControls','UseDoseControls'};
nTypes = numel(UseNames);
for ii = 1:nTypes
    
    % Replicate Use* arrays to match number of experiments
    isSeed = strcmp(UseNames{ii}, 'UseSeeds');
    if isSeed
        opts.(UseNames{ii}) = repmat(opts.(UseNames{ii})(:,1), 1, nCon);
    else
        opts.(UseNames{ii}) = repmat(opts.(UseNames{ii})(1), nCon, 1);
    end
    
    % Choose subsets of parameters to fit for each experiment
    for i_con = 1:nCon
        if isSeed
            opts.(UseNames{ii})(:,i_con) = selectUses(opts.(UseNames{ii})(:,i_con));
        else
            opts.(UseNames{ii}){i_con} = selectUses(opts.(UseNames{ii}){i_con});
        end
    end
    
end

% Generate random values for parameters by varying provided parameters up
% and down two-fold
origValues = m.k(opts.UseParams);
randomMults = 2.^(2*rand(sum(opts.UseParams),1)-1);
mOut = m;
mOut.k(opts.UseParams) = randomMults.*origValues;
mOut = mOut.Update(mOut.k);

UseTypes = {'s','q','h'};
conOut = repmat(con, nCon, 1);
for ii = 1:nTypes
    
    UseType = UseTypes{ii};
    UseName = UseNames{ii};
    isSeed = strcmp(UseName, 'UseSeeds');
    
    for i_con = nCon:-1:1
        
        % Pull out use*
        if isSeed
            Use_di = opts.(UseName)(:,i_con);
        else
            Use_di = opts.(UseName){i_con};
        end
        
        % Get original values for parameters of this type
        origValues = con.(UseType)(Use_di);
        
        randomMults = 2.^(2*rand(sum(Use_di),1)-1);
        
        if any(Use_di) % Avoids problems with empties when Use_di is all false
            conOut(i_con).(UseType)(Use_di) = randomMults.*origValues;
        end
        
        % Update experiments on last parameter type
        if ii == nTypes
            conOut(i_con) = conOut(i_con).Update(...
                conOut(i_con).s,...
                conOut(i_con).q,...
                conOut(i_con).h);
        end
        
    end
end

% Determine random time points, experiments, and outputs for the objective
% function
timelist = tF*rand(nMeasurements, 1);
outputlist = randi(m.ny, nMeasurements, 1);
experimentlist = randi(nCon, nMeasurements, 1);

sd = sdLinear(0.1, 0.01);

% Get observations
for iCon = nCon:-1:1
    isExperiment = experimentlist == iCon;
    obs(iCon) = observationLinearWeightedSumOfSquares(...
        outputlist(isExperiment),...
        timelist(isExperiment),...
        sd,...
        sprintf('Experiment %d', iCon));
end

% Simulate new experiments to get measurement values
% Random noise is added according to sd when calling sim.measurements
sim = SimulateSystem(mOut, conOut, obs, opts);
measurements = {sim.measurements};

% Create objective functions
obj = objectiveZero([nCon 1]);
for iCon = nCon:-1:1
    obj(iCon) = obs(iCon).Objective(measurements{iCon});
end

end