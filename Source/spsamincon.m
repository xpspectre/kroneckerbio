function [X,FVAL,EXITFLAG,OUTPUT,LAMBDA,GRAD,HESSIAN] = spsamincon(problem, opts)
% Optimize using simultaneous perturbation stochastic approximation algorithm.
% Exports an interface similar to the updated fmincon interface.
% Note: only supports box constraints lb and ub for now
%   Need other types of constraints
% NOTE: currently doesn't really work - need to balance 2 things:
%   small steps so variances don't go out of bounds
%   large steps so it actually converges

% Additional spsa params
if nargin < 2
    opts = [];
end

% SPSA specific params
% Guidelines for coefficients and gain sequence parameters:
%   high noise: small a, large c; opposite for low noise
%   asymptotically optimal alpha = 1.0, gamma = 1/6
%   selected values below are lowest/fastest values still satifying theory
%       large amount of data -> use larger values closer to optimal
% Note: hardcode vals in here for example/testing, make a function that applies
% guidelines automatically later. May need a couple of initial steps to estimate
% values for these guidelines (a kind of "burn in")
opts_ = [];
opts_.a = 0.1; % a/(A+1)^alpha * gkhat ~ change in magnitude of theta in early steps
opts_.A = 100; % <=10% of the max iterations
opts_.c = 2; % std dev of measurement noise, else a small positive noise floor for Bernoulli +/-1
opts_.alpha = 0.602;
opts_.gamma = 0.101;

opts = mergestruct(opts_, opts);

a     = opts.a;
A     = opts.A;
c     = opts.c;
alpha = opts.alpha;
gamma = opts.gamma;

% General/fmincon-like params
maxIter   = problem.options.MaxIter;
disp      = problem.options.Display;
outputFcn = problem.options.OutputFcn;
objLimit  = problem.options.ObjectiveLimit;

lb = problem.lb;
ub = problem.ub;

if ~isempty(problem.Aineq) || ~isempty(problem.bineq) || ~isempty(problem.Aeq) || ~isempty(problem.beq)
    warning('spsamincon:invalidConstraints', 'Non box constraints found. Ignoring.')
end

% Starting parameter vector
x = problem.x0;
nx = length(x);

% Make sure starting parameter vector falls inside bounds, or warn and fix
x = fixBounds(x, lb, ub);

% Objective function value
objective = problem.objective;

% % Print some starting summary information
y = objective(x);
fprintf('Starting objective function value: %g\n', y);

% For outputFcn
optimValues = [];
optimValues.fval = y;

% Blocking
threshold = 25;

% Keep track of obj fun vals over time
xList = x;

for k = 1:maxIter
    
    while true % Keep trying until obj fun val decreases enough/doesn't increase too much
        
        fprintf('Iter = %g\t', k);
        
        % Make random perturbation vector
        ak = a/(A + k)^alpha;
        ck = c/k^gamma;
        deltak = bernoulli(nx,1);
        
        % Make new search points
        x1 = x + ck*deltak;
        x2 = x - ck*deltak;
        
        % Make sure points fall inside bounds
        x1 = fixBounds(x1, lb, ub);
        x2 = fixBounds(x2, lb, ub);
        
        % TODO: Additional constraints
        % For NLME, make sure Sigma and Omega are valid covariance matrices, positive semi-definite
        
        % Evaluate at the 2 points
        y1 = objective(x1);
        y2 = objective(x2);
        
        % Sanity checks on values
        % Note: returns complex numbers sometimes for invalid covariance matrix
        % values in NLME, for example
        
        % Calculate gradient
        ghatk = (y1 - y2)/(2*ck);
        step = ak*ghatk;
        
        fprintf('deltaG = %g\tstep = %g\t', ghatk, ak*ghatk);
        
        % Update parameters
        x = x - step;
        
        % Show obj fun val at this step
        % Note: not necessary but useful
        y = objective(x);
        fprintf('G = %g\t', y);
        
        % See if we're done
        % optimValues.fval = objective function value
        optimValues.fval = y;
        if outputFcn(x, optimValues, [])
            return
        end
        
        % Blocking in 2nd order SPSA
        %   Discard step if new value is significantly worse than last value
        if (x - xList(end)) > threshold
            fprintf('New G too large, retrying\n');
            x = x + step; % reset this
        else
            fprintf('\n');
            xList = [xList; x];
            break
        end
            
    end
    
end

fprintf('Done.\n');

X        = x;
FVAL     = y;
EXITFLAG = 0;
OUTPUT   = [];
LAMBDA   = [];
GRAD     = [];
HESSIAN  = [];

end

function out = bernoulli(xdims,ydims)
% Generate xdims x ydims matrix of Bernoulli +/- 1 random variables
r = rand(xdims,ydims);
out = (r > 0.5) - (r <= 0.5);
end

function x = fixBounds(x, lb, ub)
% Make sure input vector falls inside box constraints
x(x<lb) = lb(x<lb);
x(x>ub) = ub(x>ub);
end

function valid = hasValidCovarianceMatrices(x,spec)
% Returns whether or not the parameter vector's params that represent elements
% of covariance matrices according to spec are valid elements of covariance
% matrices.
valid = true;
end