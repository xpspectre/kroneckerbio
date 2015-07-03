function newCell = combineCellMats(baseCell, topCell, mode)
% Combine two cell matrices of double vectors, specifying overwriting or summing
% baseCell with topCell values 
% Inputs:
%   baseCell [ 1 x nCell cell matrix of nX x 1 double vectors ]
%   topCell [ 1 x nCell cell matrix of nX x 1 double vectors ]
%   mode [ string { 'overwrite' } ]
%       Combination mode allowing:
%           'overwrite': replace values in baseCell with nonzero values in topCell
%           'sum': sum values in double vectors
% Outputs:
%   newCell [ 1 x nCell cell matrix of nX x 1 double vectors ]

if nargin < 3
    mode = 'overwrite';
end

assert(size(baseCell,1) == 1)
nCell = length(baseCell);

newCell = cell(size(baseCell));
for i = 1:nCell
    baseCell_i = baseCell{i};
    topCell_i = topCell{i};
    assert(all(size(baseCell_i) == size(topCell_i)))
    
    if strcmpi(mode, 'overwrite')
        newCell_i = baseCell_i.*logical(~topCell_i) + topCell_i.*logical(topCell_i);
    elseif strcmpi(mode, 'sum')
        newCell_i = baseCell_i + topCell_i;
    end
    
    newCell{i} = newCell_i;
end