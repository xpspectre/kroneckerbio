function [m, conNew] = updateAll(m, con, T, UseParams, UseSeeds, UseInputControls, UseDoseControls)
%updateAll Update Kronecker Bio structures when parameters change
%
%   [m, con] = updateAll(m, con, obj, T, UseParams, UseSeeds, UseInputControls, UseDoseControls)
%
%   This function performs the often needed task of updating m, con, and
%   obj whenever the active parameters T change.

% Handle vectors of models
% Below, assume all models have the same structure but different k and
%   associated derivatives
nTop = length(m);

% Constants
nk = m(1).nk;
ns = m(1).ns;
nCon = size(con, 1);
nTk = getnTk(UseParams);
nTs = nnz(UseSeeds);
[UseInputControls, nTq] = fixUseControls(UseInputControls, nCon, cat(1,con.nq));
[UseDoseControls, nTh] = fixUseControls(UseDoseControls, nCon, cat(1,con.nh));
nT = nTk + nTs + nTq + nTh;

% Update parameter sets
k = kVec2k(T(1:nTk), UseParams);

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
%   Note: currently, k is a matrix of all k's in cons
for i = 1:nTop
    kMask_i = ~isnan(k(:,i)); % mask for fit k in model i
    kFit_i = k(:,i);
    kFit_i(isnan(kFit_i)) = 0;
    ki = kFit_i.*kMask_i + m(i).k.*~kMask_i;
    m(i) = m(i).Update(ki);
end

% Update experimental conditions
if ~isnumeric(con)
    conNew = experimentZero(nCon);
%     if all(size(vec(k)) == size(k)) % expand k if there's only 1 model
%         k = repmat(k,1,nCon);
%     end
    for i = 1:nCon
        conNew(i) = con(i).Update(k(:,i), s(:,i), q{i}, h{i});
    end
end
