function F = computeObjInfo(m, con, obj, opts, dxdTSol)
% Compute Fisher information matrix of objective function w.r.t. params used in condition and
% model
%
% Inputs:
%   m [ model struct ]
%       Single model
%   con [ condition struct ]
%       Single condition built off model
%   obj [ nObj vector of objective structs ]
%       Objective functions belonging to this condition
%   opts [ options struct ]
%       Options struct
%   dxdTSol [ sensitivity solution struct ]
%       A structure vector containing the solution to the model
%       sensitivities can be provided to prevent this
%       method from recalculating it
%
% Outputs:
%   F [ 1 x 1 cell of 4 x 4 cell vector of nXi x nXj double vectors ]
%       Objective function Hessian w.r.t. params used in m and con in standard
%       Tlocal form
%
% Note: this function has been simplified to only allow a single condition

if nargin < 5
    dxdTSol = [];
end

% Process options
% Note: this currently acts as a shim to existing code - a lot of this can be
% cleaned up in the future
opts.UseParams        = logical(con.Extra.ParamsSpec{1});
opts.UseSeeds         = logical(con.Extra.ParamsSpec{2});
opts.UseInputControls = logical(con.Extra.ParamsSpec{3});
opts.UseDoseControls  = logical(con.Extra.ParamsSpec{4});

opts.continuous = any(obj(:).Continuous);
opts.RelTol = con.Extra.RelTol;
AbsTol = fixAbsTol(con.Extra.AbsTol, 3, opts.continuous, m.nx, 1, opts.UseAdjoint, opts.UseParams, opts.UseSeeds, {opts.UseInputControls}, {opts.UseDoseControls});
opts.AbsTol = AbsTol{1};
opts.ObjWeights = [obj(:).Weight];

% Constants
nx = m.nx;
ny = m.ny;
nTk = nnz(opts.UseParams);
nTs = nnz(opts.UseSeeds);
nTq = nnz(opts.UseInputControls);
nTh = nnz(opts.UseDoseControls);
nT = nTk + nTs + nTq + nTh;
paramMapper = ParamMapperOneModelType(con);
nObj = length(obj);

% Integrate
if isempty(dxdTSol)
    if opts.Verbose; fprintf('Integrating sensitivity:\n'); end
    verbose_all = max(opts.Verbose-1,0);
    if verbose_all; tic; end
    ints = integrateAllSens(m, con, obj, opts);
else
    ints = dxdTSol;
end

% Compute FIM
F = zeros(nT,nT); % T_T
for i = 1:nObj
    F = F + obj(i).F(ints(i));
end

% Print out something verbose
% if opts.Verbose; fprintf('Summary: |dGdT| = %g\t||d2GdT2|| = %g\n', norm(D), det(H)); end

% Convert FIM to standard form for local T
F = paramMapper.T2Tlocal(F);
