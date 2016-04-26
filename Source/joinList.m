function list = joinList(chunks)
% Join chunks back together into a single list. Naive implementation that just
%   concatenates chunks together w/o preallocating - this is fine when joining
%   Composite objects split on a local parpool. Note that Composite objects are
%   stored as row vectors (instead of the col vectors that splitList outputs),
%   but this is handled.
% TODO: Preallocate, considering type of things inside chunks

list = chunks{1};
nChunks = length(chunks);
for i = 2:nChunks
    list = [list; chunks{i}];
end

%% More efficient implementation
% % Get bound on the total number of elements from list
% nChunks = size(chunks,1);
% nItemsPerChunk = size(chunks{1},1);
% maxItems = nChunks * nItemsPerChunk;
% 
