function [fit, G, D] = FitObjectiveNlme(fit, opts)
%FitObjectiveNlme Optimize the Theta, Omega, and Sigma of a population nonlinear
%mixed effects (NLME) model
%
%   [fit, G, D] = FitObjectiveNLME(fit, opts)
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
%   The formulation is based on Almquist, J., Leander, J., & Jirstrand, M.
%       (2015). Using sensitivity equations for computing gradients of the
%       FOCE and FOCEI approximations to the population likelihood. Journal
%       of Pharmacokinetics and Pharmacodynamics, 191–209. doi:10.1007/s10928-015-9409-1
%
%   The NLME problem finds the optimal values of the fixed effects
%   parameters {Theta, Omega, Sigma} while accounting for random effects
%   parameters Eta. Theta are the population average values, Omega is the
%   covariance matrix of the random effects parameters, which have a
%   normal(0,Omega) distribution, and Sigma is the covariance matrix of the
%   "measurement" error terms, which have a normal(0,Rij) distribution with
%   terms in Rij including sigma.
%
%   Inputs:
%       fit: [ FitObject ]
%           Fit object that contains fitting scheme
%       opts: [ struct {} ]
%           Additional options allowing the fields specified in
%           FitObject's constructor, overriding FitObject's opts
%
%   Outputs
%       fit: [ FitObject ]
%           Fit object updated with optimized values
%       G: [ real scalar ]
%           The optimum objective function value
%       D: [ real vector nT ]
%           The objective gradient at the optimum parameter set
%
%   Notes: 
%       - Only compatible with FitObject
%       - Global optimization has been disabled for now
%       - Parameter normalization doesn't work for now
%       - Requires forward sensitivity calculation (adjoint doesn't work
%       for now) - adjoint probably won't work because lots of different
%       types of sensitivities are needed
%       - Throughout, T refers to all params, and is what's fed into obj
%       fun calculators, Theta and Eta are separate and optimized on
%
%   TODO:
%       - Make a variable to keep track of etas, for better starting
%       guesses for inner optimization in FOCEI

%% Work-up
% Clean up inputs
if nargin > 1
    fit.addOptions(opts);
end

% For convenience, copy fit object's options into this space
opts = fit.options;

% Set logging verbosity
verbosity = opts.Verbose;
assert(isnumeric(verbosity) && verbosity >= 0, 'KroneckerBio:FitObjectiveNlme:InvalidVerbosity', 'opts.Verbose must be a nonnegative number.')
opts.Verbose = max(opts.Verbose-1,0);

% Ensure Restart is a positive integer
if ~(opts.Restart >= 0)
    opts.Restart = 0;
    warning('KroneckerBio:FitObjectiveNlme', 'opts.Restart was not nonegative. It has been set to 0.')
end

if ~(opts.Restart == floor(opts.Restart))
    opts.Restart = floor(opts.Restart);
    warning('KroneckerBio:FitObjectiveNlme', 'opts.Restart was not a whole number. It has been floored.')
end

% Ensure RestartJump is a function handle
if isnumeric(opts.RestartJump)
    opts.RestartJump = @(iter,G)(opts.RestartJump);
end

%% NLME-specific options
% Enable (n+1)-th order sensitivity calculation
opts.ComputeSensPlusOne = true;

% Disable parameter normalization
if opts.Normalized
    opts.Normalized = false;
    warning('KroneckerBio:FitObjectiveNlme:NormalizationDisabled', 'FitObjectiveNlme doesn''t currently support normalized parameters. opts.Normalized set to false.')
end

% Disable adjoint sensitivity calculation
if opts.UseAdjoint
   opts.UseAdjoint = false;
   warning('KroneckerBio:FitObjectiveNlme:AdjointSensitivityDisabled', 'FitObjectiveNlme doesn''t currently support adjoint sensitivity calculation. opts.UseAdjoint set to false.')
end

% Disable complex step differentiation
if opts.ComplexStep
    opts.ComplexStep = false;
    warning('KroneckerBio:FitObjectiveNlme:ComplexStepDifferentiationDisabled', 'FitObjectiveNlme doesn''t currently support complex step differentiation. opts.ComplexStep set to false.')
end

%% Outer (theta) optimization options
% Default options
outerOpts_ = [];
outerOpts_.Algorithm               = 'active-set';
outerOpts_.UseFiniteDifferences    = false;
outerOpts_.TolOptim                = 1e-4;
outerOpts_.TolX                    = 0;
outerOpts_.MaxFunEvals             = 20000; % more fun evals/step 
outerOpts_.MaxIter                 = 1000;
outerOpts_.MaxStepSize             = 0.1; % smaller
outerOpts_.RelLineSrchBndDuration  = Inf;
outerOpts_.TolCon                  = 1e-4; % larger
outerOpts_.DerivativeCheck         = 'off'; % Turn on to debug analytic gradient calculation
outerOpts_.FinDiffType             = 'central'; % Slower but more accurate

outerOpts_ = mergestruct(outerOpts_, opts);

if outerOpts_.UseFiniteDifferences
    outerOpts_.GradObj = 'off';
else
    outerOpts_.GradObj = 'on';
end

outerOpts = optimoptions('fmincon');
outerOpts.Algorithm               = outerOpts_.Algorithm;
outerOpts.TolFun                  = outerOpts_.TolOptim;
outerOpts.TolX                    = outerOpts_.TolX;
outerOpts.OutputFcn               = @isTerminalObjOuter;
outerOpts.GradObj                 = outerOpts_.GradObj;
outerOpts.Hessian                 = 'off';
outerOpts.MaxFunEvals             = outerOpts_.MaxFunEvals;
outerOpts.MaxIter                 = outerOpts_.MaxIter;
outerOpts.RelLineSrchBnd          = outerOpts_.MaxStepSize;
outerOpts.RelLineSrchBndDuration  = outerOpts_.RelLineSrchBndDuration;
outerOpts.TolCon                  = outerOpts_.TolCon;
outerOpts.DerivativeCheck         = outerOpts_.DerivativeCheck;
outerOpts.FinDiffType             = outerOpts_.FinDiffType;

if verbosity > 1
    outerOpts.Display = 'iter'; % or 'iter-detailed'
else
    outerOpts.Display = 'off';
end

%% Rewrite FitObject's options
fit.addOptions(opts);

%% Construct starting variable parameter set of thetas
[thetas, etas] = getParameterMask({fit.Models(1).Parameters.Name});
T0 = fit.collectParams;
nT = length(T0);
Theta0 = T0(thetas);
nTheta = length(Theta0);

% Collect bounds
[LowerBound, UpperBound] = fit.collectBounds;
LowerBound = LowerBound(thetas);
UpperBound = UpperBound(thetas);

% Apply bounds to starting parameters before optimizing
% fmincon will choose a weird value if a starting parameter is outside the bounds
Theta0(Theta0 < LowerBound) = LowerBound(Theta0 < LowerBound);
Theta0(Theta0 > UpperBound) = UpperBound(Theta0 > UpperBound);

% Whether to use conditional expectation
%   Enables inner optimization to find mode of eta
if strcmp(opts.Approximation, 'FOCEI')
    useCE = true;
else
    useCE = false;
end

% Number of participants
ni = fit.nConditions;

%% Abort if no fit params specified
if numel(Theta0) == 0
    [G, D] = outerObjective(Theta0);
    return
end

%% Run optimization
% Initialize loop
ThetaHat  = Theta0;
Gbest = Inf;
ThetaBest = Theta0;

for iRestart = 1:opts.Restart+1
    % Init abort parameters
    aborted = false;
    ThetaAbort = ThetaHat;
    
    % Create local optimization problem
    %   Always needed - used as a subset/refinement of global optimization
    localProblem = createOptimProblem('fmincon', 'objective', @outerObjective, ...
        'x0', ThetaHat, 'Aeq', opts.Aeq, 'beq', opts.beq, ...
        'lb', LowerBound, 'ub', UpperBound, 'options', outerOpts);
    
    if opts.Verbose; fprintf('Beginning gradient descent...\n'); end
    [ThetaHat, G, exitflag, ~, ~, D] = fmincon(localProblem);
    
    % Check abortion status
    % Abortion values are not returned by fmincon and must be retrieved
    if aborted
        ThetaHat = ThetaAbort;
        [G, D] = outerObjective(ThetaHat);
    end
    
    % Re-apply stiff bounds
    ThetaHat(ThetaHat < LowerBound) = LowerBound(ThetaHat < LowerBound);
    ThetaHat(ThetaHat > UpperBound) = UpperBound(ThetaHat > UpperBound);
    
    % See if these parameters are better than previous iterations
    if G < Gbest
        Gbest = G;
        ThetaBest = ThetaHat;
    end

    % Terminate if goal has been met
    if G <= opts.TerminalObj
        break
    end
    
    % Jump parameters before restarting
    % Retain the deterministic nature of fitting by fixing the stream
    rng_state = rng;
    ThetaHat = inf2big(ThetaHat); % fix -Infs (can occur if normalized and got 0)
    prodThat = prod(ThetaHat); % model fingerprint
    rng(mod(prodThat/eps(prodThat)*iRestart,2^32));
    ThetaHat = exp(log(ThetaBest) + randn(nTheta,1) .* vec(opts.RestartJump(iRestart,G)));
    rng(rng_state);
    
    % Prevent jumps from leaving bounds
    while any(ThetaHat < LowerBound) || any(ThetaHat > UpperBound)
        ThetaHat(ThetaHat < LowerBound) = 2*(LowerBound(ThetaHat < LowerBound)) - ThetaHat(ThetaHat < LowerBound);
        ThetaHat(ThetaHat > UpperBound) = 2*(UpperBound(ThetaHat > UpperBound)) - ThetaHat(ThetaHat > UpperBound);
    end
end

% Update parameter sets
T = zeros(nT,1);
T(thetas) = ThetaBest;
fit.updateParams(T);

% End of function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outer (theta) objective function %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [G, D] = outerObjective(Theta)
        % Objective function optimizes on vector or thetas.
        % Note: Thetas are put into the total vector of all params for
        %   computeObjective. Gradient is return w.r.t. thetas only.
        
        % If fmincon chooses values outside the bounds, force them to
        % equal to the bound
        Theta(Theta < LowerBound) = LowerBound(Theta < LowerBound);
        Theta(Theta > UpperBound) = UpperBound(Theta > UpperBound);
        
        %% Update parameter sets
        T = zeros(nT,1);
        T(thetas) = Theta;
        
        % Perform inner optimization for FOCEI
        % Calculate obj fun val and grad using optimal etas for each
        %   individual
        if useCE % FOCEI
            fprintf('Performing inner optimization...\n') % DEBUG
            
            if nargout == 1
                G = 0;
                for i = 1:ni
                    Gi = calcIndividualContribution(i);
                    G = G + Gi;
                end
            end
            if nargout == 2
                G = 0;
                D = zeros(nTheta,1);
                for i = 1:ni
                    [Gi, Di] = calcIndividualContribution(i);
                    G = G + Gi;
                    D = D + Di;
                end
            end
            
        else % CO
            fit.updateParams(T);
            
            % Compute objective (and gradient if desired)
            if nargout == 1
                G = fit.computeObjective;
            end
            if nargout == 2
                [G, D] = fit.computeObjective;
                D = D(thetas);
            end
        end
        
    end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Halt outer optimization on terminal goal %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function stop = isTerminalObjOuter(x, optimValues, state)
        if optimValues.fval < opts.TerminalObj
            aborted = true;
            ThetaAbort = x;
            stop = true;
        else
            stop = false;
        end
    end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate individual contribution in FOCEI %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [G, D] = calcIndividualContribution(i)
        
        % Create inner fitting problem
        innerFit = FitObject(['Fit_FOCEI_Inner_', num2str(i)]);
        innerFit.addModel(fit.Models(1));
        innerObs = observationIndividualVariability(fit.Objectives(i).private.outputs, fit.Objectives(i).private.times);
        innerObj = innerObs.Objective(fit.Objectives(i).private.measurements);
        innerFit.addFitConditionData(innerObj, fit.Conditions(i));
        
        % Run inner fitting problem
        [~, ~, ~, eta_i] = fitIndividualVariability(innerFit, opts);
        
        % Calculate obj fun val and grad cont
        T_i = T;
        T_i(etas) = eta_i;
        fit.updateParams(T_i); % temporarily update fit to the individual's optimal eta
        
        if nargout == 1
            G = fit.computeObjective(i);
        end
        if nargout == 2
            [G, D] = fit.computeObjective(i);
            D = D(thetas);
        end
        
        % Reset fit's params - is this inefficient?
        fit.updateParams(T);
        
    end
end
