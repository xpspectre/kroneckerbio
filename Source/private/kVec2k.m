function k = kVec2k(kVec, UseParams, conIdx)
% inverse of k2kVec

%   conIdx [ optional integer ]

% Values that aren't fit are filled with NaNs; calling function needs to merge this with nonfit k when updating

[nk, nCon] = size(UseParams);

k = NaN(nk, nCon);
startIdx = 1;
for i = 1:nk
    
    % Get appropriate positions in k
    kUse = UseParams(i,:);
    kPos = logical(kUse);
    
    % Get appropriate variable parameters
    kUse = UseParams(i,:);
    kUse(kUse==0) = [];
    nUniqueNonZero = length(unique(kUse));
    kElements = kVec(startIdx:startIdx+nUniqueNonZero-1);
    
    % Get appropriate multiplicity in k
    kUse = UseParams(i,:);
    kUse(kUse==0) = []; % NaNs are unique
    [~,~,ic] = unique(kUse, 'stable');
    
    % Assign
    k(i,kPos) = kElements(ic);
    
    % Increment counter for position in kVec
    startIdx = startIdx + nUniqueNonZero;
    
end

if nargin > 2
    k = k(:,conIdx);
end
