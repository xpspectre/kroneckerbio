function G = computeObj(m, con, obj, opts)
% Compute objective function value w.r.t. params used in condition and model
% Inputs:
%   m [ model struct ]
%       Single model
%   con [ condition struct ]
%       Single condition built off model
%   obj [ nObj vector of objective structs ]
%       Objective functions belonging to this condition
%   opts [ struct ]
%       Options struct
%       
% Outputs:
%   G [ scalar double ]
%       Objective function value at params specified in m and con
%
% Note: this function has been simplified to only allow a single condition
verbose = logical(opts.Verbose);
verbose_all = max(verbose-1,0);

% Process options
% Note: this currently acts as a shim to existing code - a lot of this can be
% cleaned up in the future
opts.UseParams        = logical(con.ParamsSpec{1});
opts.UseSeeds         = logical(con.ParamsSpec{2});
opts.UseInputControls = logical(con.ParamsSpec{3});
opts.UseDoseControls  = logical(con.ParamsSpec{4});

opts.continuous = any(obj(:).Continuous);
opts.RelTol = con.RelTol;
AbsTol = fixAbsTol(con.AbsTol, 1, opts.continuous, m.nx, 1, opts.UseAdjoint, opts.UseParams, opts.UseSeeds, {opts.UseInputControls}, {opts.UseDoseControls});
opts.AbsTol = AbsTol{1};
opts.ObjWeights = [obj(:).Weight];

% Constants
nObj = length(obj);

% Initialize variables
G = 0;

if verbose; disp('Integrating forward...'); end
if verbose_all; tic; end

% Integrate
ints = integrateAllSys(m, con, obj, opts);

% Add fields for prior objectives
ints.UseParams        = opts.UseParams;
ints.UseSeeds         = opts.UseSeeds;
ints.UseInputControls = opts.UseInputControls;
ints.UseDoseControls  = opts.UseDoseControls;

% Extract continuous term
if opts.continuous
    error('Continuous objective functions have not been updated')
    nx = m.nx;
    G_cont = ints.y(nx+1,end);
else
    G_cont = 0;
end

% Compute discrete term
G_disc = 0;
for i = 1:nObj
    G_disc_i = obj(i).G(ints(i));
    G_disc = G_disc + opts.ObjWeights(i) * G_disc_i;
end

% Add to cumulative goal value
G = G + G_cont + G_disc;

if verbose_all; fprintf('G = %g\tTime = %0.2f\n', G_cont + G_disc, toc); end

if verbose; fprintf('Summary: G = %g\n', G); end
