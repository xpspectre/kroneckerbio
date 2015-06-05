function kMap = getkMap(UseParams)

% Gets position of components in full k matrix in vector of fit params kT

% Make map of position in k to kVec
%   0 indicates not fit
[nk, nCon] = size(UseParams);
kMap = zeros(nk, nCon);
rowStart = 0;
for i = 1:nk
    row = UseParams(i,:);
    row(row==0) = [];
    [c, ~, ic] = unique(row, 'stable');
    rowPos = rowStart + ic;
    kMap(i,UseParams(i,:)~=0) = rowPos;
    rowStart = rowStart + length(c);
end
