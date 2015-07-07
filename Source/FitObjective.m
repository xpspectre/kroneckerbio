function varargout = FitObjective(varargin)
%FitObjective Optimize the parameters of a model to minimize a set of
%   objective functions
%
%   Mathematically: T = argmin(G(T))
%
%   [fit, G, D] = FitObjective(fit)
%   [fit, G, D] = FitObjective(fit, opts)
%   [m, con, G, D] = FitObjective(m, con, obj)
%   [m, con, G, D] = FitObjective(m, con, obj, opts)
%
%   FitObjective uses the derivatives in the Kronecker model and in the
%   objective functions to build a function that can not only evaluate the
%   objective functions at particular parameter sets, but also evaluate the
%   gradient at those parameter sets, thereby pointing in the direction of
%   a more optimum parameter set. This function is built around Matlab's
%   fmincon, which is a gradient descent minimizer. It varies the
%   parameters attempting to find the parameter set that will minimize the
%   objective function while keeping with the bounds.
%
%   Inputs:
%       fit: [ FitObject ]
%           Fit object that contains fitting scheme
%       m: [ model struct scalar ]
%           The KroneckerBio model that will be simulated
%       con: [ experiment struct vector ]
%           The experimental conditions under which the model will be simulated
%       obj: [ objective struct matrix ]
%           The objective structures defining the objective functions to be
%           evaluated.
%       opts: [ struct {} ]
%           Options struct allowing the fields specified in FitObject's constructor
%
%   Outputs
%       fit: [ FitObject ]
%           Fit object updated with optimized values
%       m: [ model scalar ]
%           The model with optimal rate parameters applied
%       con: [ experiment vector ]
%           The experimental conditions with optimal seed, input control, and
%           dose control parameters applied
%       G: [ real scalar ]
%           The optimum objective function value
%       D: [ real vector nT ]
%           The objective gradient at the optimum parameter set
%
%   Generally compatible with old (single model, vector of conditions, matrix of
%   objectives) and new (FitObject) forms of specifying fitting scheme. Outputs
%   results in the same format as input.

% (c) 2015 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

%% Work-up
% Clean up inputs
switch nargin
    case 1
        fit = varargin{1};
    case 2
        fit = varargin{1};
        fit.addOptions(varargin{2});
    case 3
        fit = FitObject.buildFitObject(varargin{1}, varargin{2}, varargin{3});
    case 4
        fit = FitObject.buildFitObject(varargin{1}, varargin{2}, varargin{3}, varargin{4});
    otherwise
        
end

% For convenience, copy fit object's options into this space
opts = fit.options;

verbose = logical(opts.Verbose);
opts.Verbose = max(opts.Verbose-1,0);

% Ensure Restart is a positive integer
if ~(opts.Restart >= 0)
    opts.Restart = 0;
    warning('KroneckerBio:FitObjective', 'opts.Restart was not nonegative. It has been set to 0.')
end

if ~(opts.Restart == floor(opts.Restart))
    opts.Restart = floor(opts.Restart);
    warning('KroneckerBio:FitObjective', 'opts.Restart was not a whole number. It has been floored.')
end

% Ensure RestartJump is a function handle
if isnumeric(opts.RestartJump)
    opts.RestartJump = @(iter,G)(opts.RestartJump);
end

% Collect bounds
[LowerBound, UpperBound] = fit.collectBounds;

% Construct starting variable parameter set
T0 = fit.collectParams;
nT = length(T0);

%% Local optimization options
localOpts = optimoptions('fmincon');
localOpts.Algorithm               = opts.Algorithm;
localOpts.TolFun                  = opts.TolOptim;
localOpts.TolX                    = 0;
localOpts.OutputFcn               = @isTerminalObj;
localOpts.GradObj                 = 'on';
localOpts.Hessian                 = 'off'; % unused
localOpts.MaxFunEvals             = opts.MaxFunEvals;
localOpts.MaxIter                 = opts.MaxIter;
localOpts.RelLineSrchBnd          = opts.MaxStepSize;
localOpts.RelLineSrchBndDuration  = Inf;
localOpts.TolCon                  = 1e-6;

if verbose
    localOpts.Display = 'iter';
else
    localOpts.Display = 'off';
end

%% Global optimization options
% TODO: make sure options are relevant for solver; extend global solver API
globalOpts = opts.GlobalOpts;

% Sanity checking
if strcmpi(globalOpts.Algorithm, 'multistart') && globalOpts.UseParallel
    warning('KroneckerBio:FitObjective', 'Using multistart with UseParallel is not supported at this time (due to global variable in obj fun usage).')
end

%% Normalize parameters
if opts.Normalized
    % Normalize starting parameters and bounds
    T0 = log(T0);
    LowerBound = log(LowerBound);
    UpperBound = log(UpperBound);
    
    % Change relative line search bound to an absolute scale in log space
    % Because fmincon lacks an absolute option, this hack circumvents that
    localOpts.TypicalX = zeros(nT,1) + log(1 + opts.MaxStepSize)*log(realmax);
    localOpts.RelLineSrchBnd = 1 / log(realmax);
end

%% Apply bounds to starting parameters before optimizing
% fmincon will choose a weird value if a starting parameter is outside the bounds
T0(T0 < LowerBound) = LowerBound(T0 < LowerBound);
T0(T0 > UpperBound) = UpperBound(T0 > UpperBound);

%% Abort in rare case of no optimization
if numel(T0) == 0
    [G, D] = objective(T0);
    return
end

%% Run optimization
% Initialize loop
That = T0;
Gbest = Inf;
Tbest = T0;

for iRestart = 1:opts.Restart+1
    % Init abort parameters
    aborted = false;
    Tabort = That;
    
    % Create local optimization problem
    %   Always needed - used as a subset/refinement of global optimization
    localProblem = createOptimProblem('fmincon', 'objective', @objective, ...
        'x0', That, 'Aeq', opts.Aeq, 'beq', opts.beq, ...
        'lb', LowerBound, 'ub', UpperBound, 'options', localOpts);
    
    % Run specified optimization
    if opts.GlobalOptimization
        
        if opts.Verbose
            fprintf('Beginning global optimization with %s...\n', globalOpts.Algorithm)
        end
        
        switch globalOpts.Algorithm
            case 'globalsearch'
                gs = GlobalSearch('StartPointsToRun', globalOpts.StartPointsToRun);
                [That, G, exitflag] = run(gs, localProblem);
            case 'multistart'
                ms = MultiStart('StartPointsToRun', globalOpts.StartPointsToRun, ...
                    'UseParallel', globalOpts.UseParallel);
                [That, G, exitflag] = run(ms, localProblem, globalOpts.nStartPoints);
            case 'patternsearch'
                psOpts = psoptimset('MaxIter', globalOpts.MaxIter, 'UseParallel', globalOpts.UseParallel);
                [That, G, exitflag] = patternsearch(@objective, That, [], [], ...
                    opts.Aeq, opts.beq, LowerBound, UpperBound, [], psOpts);
            otherwise
                error('Error:KroneckerBio:FitObjective: %s global optimization algorithm not recognized.', globalOpts.Algorithm)
        end
        
        [~, D] = objective(That); % since global solvers don't return gradient at endpoint
        
    else
        
        if opts.Verbose; fprintf('Beginning gradient descent...\n'); end
        [That, G, exitflag, ~, ~, D] = fmincon(localProblem);
        
    end
    
    % Check abortion status
    % Abortion values are not returned by fmincon and must be retrieved
    if aborted
        That = Tabort;
        [G, D] = objective(That);
    end
    
    % Re-apply stiff bounds
    That(That < LowerBound) = LowerBound(That < LowerBound);
    That(That > UpperBound) = UpperBound(That > UpperBound);
    
    % See if these parameters are better than previous iterations
    if G < Gbest
        Gbest = G;
        Tbest = That;
    end

    % Terminate if goal has been met
    if G <= opts.TerminalObj
        break
    end
    
    % Jump parameters before restarting
    % Retain the deterministic nature of fitting by fixing the stream
    rng_state = rng;
    That = inf2big(That); % fix -Infs (can occur if normalized and got 0)
    prodThat = prod(That); % model fingerprint
    rng(mod(prodThat/eps(prodThat)*iRestart,2^32));
    if opts.Normalized
        That = Tbest + randn(nT,1) .* vec(opts.RestartJump(iRestart,G));
    else
        That = exp(log(Tbest) + randn(nT,1) .* vec(opts.RestartJump(iRestart,G)));
    end
    rng(rng_state);
    
    % Prevent jumps from leaving bounds
    while any(That < LowerBound) || any(That > UpperBound)
        That(That < LowerBound) = 2*(LowerBound(That < LowerBound)) - That(That < LowerBound);
        That(That > UpperBound) = 2*(UpperBound(That > UpperBound)) - That(That > UpperBound);
    end
end

% Unnormalize
if opts.Normalized
    Tbest = exp(Tbest);
    D = D ./ Tbest;
end

% Update parameter sets
fit.updateParams(Tbest);

% Output according to input signature
switch nargin
    case 1
        varargout{1} = fit;
        varargout{2} = G;
        varargout{3} = D;
    case 2
        varargout{1} = fit;
        varargout{2} = G;
        varargout{3} = D;
    case 3
        varargout{1} = fit.Models;
        varargout{2} = fit.Conditions;
        varargout{3} = G;
        varargout{4} = D;
    case 4
        varargout{1} = fit.Models;
        varargout{2} = fit.Conditions;
        varargout{3} = G;
        varargout{4} = D;
    otherwise
        
end

% End of function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Objective function for fmincon %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [G, D] = objective(T)
        
        % Unnormalize
        if opts.Normalized
            T = exp(T);
        else
            % If fmincon chooses values that violate the lower bounds, force them to be equal to the lower bounds
            T(T < LowerBound) = LowerBound(T < LowerBound);
            T(T > UpperBound) = UpperBound(T > UpperBound);
        end
        
        % Update parameter sets
        fit.updateParams(T);
        
        if nargout == 1
            G = fit.computeObjective;
        end
        if nargout == 2
            [G, D] = fit.computeObjective;
            
            % Normalize gradient
            if opts.Normalized
                D = D .* fit.collectParams;
            end
        end
        
    end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Halt optimization on terminal goal %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function stop = isTerminalObj(x, optimValues, state)
        if optimValues.fval  < opts.TerminalObj
            aborted = true;
            Tabort = x;
            stop = true;
        else
            stop = false;
        end
    end

end

