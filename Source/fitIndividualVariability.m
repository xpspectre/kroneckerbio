function [fit, G, D, EtaBest] = fitIndividualVariability(fit, opts)
%fitIndividualVariabilities Optimize the etas of an individual

%% Work-up
% Clean up inputs
if nargin > 1
    fit.addOptions(opts);
end

% For convenience, copy fit object's options into this space
opts = fit.options;

% Set logging verbosity
verbosity = opts.Verbose;
assert(isnumeric(verbosity) && verbosity >= 0, 'KroneckerBio:fitIndividualVariability:InvalidVerbosity', 'opts.Verbose must be a nonnegative number.')
opts.Verbose = max(opts.Verbose-1,0);

%% NLME-specific options
% Enable (n+1)-th order sensitivity calculation
% NOTE/TODO: eta obj fun val and grad both require only up to 1st order sensitivities)
opts.ComputeSensPlusOne = true;

% Disable parameter normalization
if opts.Normalized
    opts.Normalized = false;
    warning('KroneckerBio:fitIndividualVariability:NormalizationDisabled', 'FitObjectiveNlme doesn''t currently support normalized parameters. opts.Normalized set to false.')
end

% Disable adjoint sensitivity calculation
if opts.UseAdjoint
   opts.UseAdjoint = false;
   warning('KroneckerBio:fitIndividualVariability:AdjointSensitivityDisabled', 'FitObjectiveNlme doesn''t currently support adjoint sensitivity calculation. opts.UseAdjoint set to false.')
end

% Disable complex step differentiation
if opts.ComplexStep
    opts.ComplexStep = false;
    warning('KroneckerBio:fitIndividualVariability:ComplexStepDifferentiationDisabled', 'FitObjectiveNlme doesn''t currently support complex step differentiation. opts.ComplexStep set to false.')
end

%% Inner (eta) optimization options
% Default options
innerOpts_ = [];
% innerOpts_.Algorithm           = 'trust-region';
innerOpts_.Algorithm             = 'quasi-newton'; % Commonly used in other implementations
innerOpts_.UseFiniteDifferences  = false;
innerOpts_.TolOptim              = 1e-5;
innerOpts_.TolX                  = 0;
innerOpts_.MaxFunEvals           = 500;
innerOpts_.MaxIter               = 100;
innerOpts_.DerivativeCheck       = 'off'; % Turn on to debug analytic gradient calculation

innerOpts_ = mergestruct(innerOpts_, opts.InnerOpts);

if innerOpts_.UseFiniteDifferences
    innerOpts_.GradObj = 'off';
else
    innerOpts_.GradObj = 'on';
end

innerOpts = optimoptions('fminunc');
innerOpts.Algorithm              = innerOpts_.Algorithm;
innerOpts.TolFun                 = innerOpts_.TolOptim;
innerOpts.TolX                   = innerOpts_.TolX;
innerOpts.OutputFcn              = @isTerminalObjInner;
innerOpts.GradObj                = innerOpts_.GradObj;
innerOpts.Hessian                = 'off';
innerOpts.MaxFunEvals            = innerOpts_.MaxFunEvals;
innerOpts.MaxIter                = innerOpts_.MaxIter;
innerOpts.DerivativeCheck        = innerOpts_.DerivativeCheck;

if verbosity > 1
    innerOpts.Display = 'iter'; % or 'iter-detailed'
else
    innerOpts.Display = 'off';
end

%% Rewrite FitObject's options
fit.addOptions(opts);

%% Construct starting variable parameter set of etas
[~, etas] = getParameterMask({fit.Models(1).Parameters.Name});
T0 = fit.collectParams;
nT = length(T0);
Eta0 = T0(etas);
nEta = length(Eta0);

%% Abort if no fit params specified
if nEta == 0
    [G, D] = innerObjective(Eta0);
    return
end

%% Run optimization
EtaHat  = Eta0;
innerProblem = createOptimProblem('fminunc', 'objective', @innerObjective, ...
    'x0', EtaHat, 'options', innerOpts);

if opts.Verbose; fprintf('Beginning gradient descent...\n'); end
[EtaBest, G, exitflag, ~, ~, D] = fminunc(innerProblem);

%% Output final results
if nargout < 4
    T = T0;
    T(etas) = EtaBest;
    fit.updateParams(T);
elseif nargout == 4
    % pass, just returning EtaBest
end

% % DEBUG - display optimal eta
% fprintf('Optimal eta is %g\n', EtaBest)
% % DEBUG - calculate obj fun val for several etas
% ntests = 100;
% testEtas = linspace(-1,1,ntests);
% vals = zeros(size(testEtas));
% for i = 1:ntests
%     vals(i) = innerObjective(testEtas(i));
% end
% figure
% plot(testEtas, vals)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inner (eta) objective function %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [G, D] = innerObjective(Eta)
        % Objective function optimizes on vector or etas.
        % Note: Etas are put into the total vector of all params for
        %   computeObjective. Gradient is return w.r.t. etas only.
        
        %% Update parameter sets
        T = T0;
        T(etas) = Eta;
        
        fit.updateParams(T);
        
        %% Compute objective (and gradient if desired)
        if nargout == 1
            G = fit.computeObjective;
        end
        if nargout == 2
            [G, D] = fit.computeObjective;
            D = D(etas);
        end
        
    end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Halt inner optimization on terminal goal %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function stop = isTerminalObjInner(x, optimValues, state)
        if optimValues.fval < opts.TerminalObj
            aborted = true;
            ThetaAbort = x;
            stop = true;
        else
            stop = false;
        end
    end



end

