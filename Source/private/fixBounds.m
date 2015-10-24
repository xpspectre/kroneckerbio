function bounds = fixBounds(bounds, UseParams, UseSeeds, UseInputControls, UseDoseControls)
% Standardize bounds to match vector of params to be optimized as expected by fmincon.
% Possible input bounds formats:
%   nT : all the fit parameters in the right order
%   nk*nCon+ns*nCon+nq+nh : all params, including those not fit
%   nk*nCon+ns*nCon : all rate and seed params, including those not fit
%   nk : rate params only, including those not fit, 1 set for all conditions
%   1 : for everything
% Need to be careful on number of bounds
% Note that parameters are assumed to be expanded out to multiple conditions before this

% Get indices of varying parameters
[~, paramMask] = k2kVec(UseParams, UseParams);

% Constants
nk = size(UseParams, 1);
ns = size(UseSeeds, 1);
nq = numel(cat(1,UseInputControls{:}));
nh = numel(cat(1,UseDoseControls{:}));
nTk = nnz(paramMask);
nTs = nnz(UseSeeds);
nTq = nnz(cat(1,UseInputControls{:}));
nTh = nnz(cat(1,UseDoseControls{:}));
nT = nTk + nTs + nTq + nTh;
nCon = size(UseSeeds, 2);

l = numel(bounds);

if l == 1
    bounds = zeros(nT,1) + bounds;
elseif l == nk*nCon+ns*nCon+nq+nh
    bounds = bounds([paramMask; vec(UseSeeds); cat(1,UseInputControls{:}); cat(1,UseDoseControls{:})]);
elseif l == nT
    %bounds = bounds; % do nothing
elseif l == nk && nTs == 0 && nTq == 0 && nTh == 0
    bounds = bounds(logical(UseParams(:,1)));
elseif l == nk*nCon+ns*nCon && nTq == 0 && nTh == 0
    bounds = bounds([paramMask; vec(UseSeeds)]);
else
    error('KroneckerBio:BoundSize', ...
        'LowerBound and UpperBound must be vectors the length of Theta, number of fit parameters; m.nk*nCon+m.ns*nCon+m.nq+m.nh, number of total parameters; m.nk*nCon+m.ns*nCon if there are no inputs or doses; m.nk if only fitting rate parameters; or scalar')
end
