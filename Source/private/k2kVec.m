function [kVec, uniqueMask] = k2kVec(k, UseParams)
% Returns vector of unique k's from input matrix of k's according to a UseParams
%   specification matrix. If multiple entries in k are different but forced to
%   be the same according to UseParams, use the 1st entry (leftmost)
% Input:
%   k [ double matrix nk by nCon ]
%   UseParams [ integer matrix nk by nCon ]
% Output:
%   kvec [ double vector nUniquek ]
%       Format: col vector of variants of k among experiemnts, then next k
%   uniqueMask [ logical vector nk*nCon ]

[nk, nCon] = size(UseParams);

kMask = NaN(nk, nCon);
for i = 1:nk
    
    % Get unique positions of k
    kUse = UseParams(i,:);
    kUse(kUse==0) = NaN; % NaNs are unique
    [~,ia,~] = unique(kUse, 'stable');
    kUse = kUse(ia)';
    kUse(isnan(kUse)) = 0;
    kPos = logical(kUse);
    
    kCons = k(i,:);
    kCons = kCons(kPos);
    
    kMask(i,kPos) = kCons;
    
end
kMask = reshape(kMask', nk*nCon, 1);
kVec = kMask;
kVec(isnan(kVec)) = [];

if nargout > 1
    uniqueMask = ~isnan(kMask);
end