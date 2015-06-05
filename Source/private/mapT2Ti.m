function T = mapT2Ti(Ti, i, opts)
% Map individual condition Ti vector to overall fmincon T vector
% Needs overall opts struct

% Sequential k refer to different conditions
% Sequential everything else is per condition

% k
UseParams_i = logical(opts.UseParams(:,i));
nki = sum(UseParams_i);
% Make dummy k with all 0's + condition i's values
k = zeros(size(opts.UseParams,1),1);
k(UseParams_i) = Ti(1:nki);
% Write elements of k matrix to kT
kMap = getkMap(opts.UseParams);
nkT = max(max(kMap));
kMap = kMap(:,i);
ki = zeros(nkT,1);
for j = 1:length(k)
    if kMap(j) ~= 0
        ki(kMap(j)) = k(j);
    end
end

% s
sStart = nnz(opts.UseSeeds(:,1:i-1))+1; % position to start at in seeds part of T
sEnd = sStart + nnz(opts.UseSeeds(:,i))-1; % last position of seeds contributed by condition i
nsi = sEnd - sStart + 1;
si = zeros(nnz(opts.UseSeeds),1);
si(sStart:sEnd) = Ti(nki+1:nki+nsi);

% q
qStart = nnz(cat(1,opts.UseInputControls{1:i-1}))+1;
qEnd = qStart + nnz(opts.UseInputControls{i})-1;
nqi = qEnd - qStart + 1;
if nqi == 0
    qi = [];
else
    qi = zeros(nnz(cat(1,opts.UseInputControls{:})),1);
    qi(qStart:qEnd) = Ti(nki+nsi+1:ki+nsi+nqi);
end

% h
hStart = nnz(cat(1,opts.UseDoseControls{1:i-1}))+1;
hEnd = hStart + nnz(opts.UseDoseControls{i})-1;
nhi = hEnd - hStart + 1;
if nhi == 0
    hi = [];
else
    hi = zeros(nnz(cat(1,opts.UseDoseControls{:})),1);
    hi(hStart:hEnd) = Ti(nki+nsi+nqi+1:ki+nsi+nqi+nhi);
end

T = [ki; si; qi; hi];