function chunks = splitList(list, n)
% Split a matrix/col vector into n ~evenly sized chunks by row, where left over elements
%   are added to first length(list) mod n chunks.
% Useful for splitting up a parallel job.
%
% Inputs:
%   list [ mx x my matrix ]
%       Matrix (usually a mx x 1 vector) containing rows to be distributed. Can
%       be any type, including a cell matrix.
%   n [ scalar positive integer ]
%       Number of chunks to distribute list into
%
% Outputs:
%   chunks [ n x 1 cell array ]
%       Split list, with each cell containing the rows of the matrix ~evenly
%       divided

% Get base chunk length
nElements = length(list);
baseChunkLength = floor(nElements/n);
remainder = mod(nElements,n);

% Make array of actual chunk lengths
lengths = zeros(1,n);
for i = 1:n
    if i <= remainder
        lengths(i) = baseChunkLength + 1;
    else
        lengths(i) = baseChunkLength;
    end
end

chunks = mat2cell(list,lengths);