function T = collectActiveParameters(m, con, UseParams, UseSeeds, UseInputControls, UseDoseControls)
% Get vector of parameters optimizer will optimize
% Gets starting k's from conditions

% Constants
nCon = size(con, 1);
nk   = m.nk;
ns   = m.ns;

% Store complete parameter sets
k = zeros(nk, nCon);
for i = 1:nCon
    k(:,i) = con(i).k;
end

s = zeros(ns, nCon);
for i = 1:nCon
    s(:,i) = con(i).s;
end

nq = sum(cat(1,con.nq));
q = zeros(nq,1);
index = 0;
for i = 1:nCon
    q(index+1:index+con(i).nq) = con(i).q;
end

nh = sum(cat(1,con.nh));
h = zeros(nh,1);
index = 0;
for i = 1:nCon
    h(index+1:index+con(i).nh) = con(i).h;
end

% Get parameters specified to be unique from each set

% Rate parameters come from m by default (when experiments are made) but can be
%   modified by experiments (TODO: unify parameters)
% If multiple experiments specify a parameter that's shared, take the 1st one
%   TODO: warn if this happens

% Construct starting variable parameter set
T = [k2kVec(k, UseParams);
     vec(s(UseSeeds));
     q(cat(1, UseInputControls{:}));
     h(cat(1, UseDoseControls{:}))];
