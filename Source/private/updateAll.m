function [m, conNew] = updateAll(m, con, T, UseParams, UseSeeds, UseInputControls, UseDoseControls)
%updateAll Update Kronecker Bio structures when parameters change
%
%   [m, con] = updateAll(m, con, obj, T, UseParams, UseSeeds, UseInputControls, UseDoseControls)
%
%   This function performs the often needed task of updating m, con, and
%   obj whenever the active parameters T change.

% Constants
nk = m.nk;
ns = m.ns;
nCon = size(con, 1);
nTk = nnz(UseParams);
nTs = nnz(UseSeeds);
[UseInputControls, nTq] = fixUseControls(UseInputControls, nCon, cat(1,con.nq));
[UseDoseControls, nTh] = fixUseControls(UseDoseControls, nCon, cat(1,con.nh));
nT = nTk + nTs + nTq + nTh;

% Update parameter sets
k = zeros(nk,nCon);
for i = 1:nCon
    k(:,i) = con(i).k;
end
k(UseParams) = T(1:nTk);

s = zeros(ns,nCon);
for i = 1:nCon
    s(:,i) = con(i).s;
end
s(UseSeeds) = T(nTk+1:nTk+nTs);

q = cell(nCon,1);
endIndex = 0;
for i = 1:nCon
    q{i} = con(i).q;
    startIndex = endIndex + 1;
    endIndex = endIndex + nnz(UseInputControls{i});
    q{i}(UseInputControls{i}) = T(nTk+nTs+startIndex:nTk+nTs+endIndex);
end

h = cell(nCon,1);
endIndex = 0;
for i = 1:nCon
    h{i} = con(i).h;
    startIndex = endIndex + 1;
    endIndex = endIndex + nnz(UseDoseControls{i});
    h{i}(UseDoseControls{i}) = T(nTk+nTs+nTq+startIndex:nTk+nTs+nTq+endIndex);
end

% Update model
%   Closures in m read this k
m = m.Update(k);

% Update experimental conditions
if ~isnumeric(con)
    conNew = experimentZero(nCon);
    for i = 1:nCon
        conNew(i) = con(i).Update(k(:,i), s(:,i), q{i}, h{i});
    end
end
