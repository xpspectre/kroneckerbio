classdef IndexAccumulator < handle
    % State machine that generates and returns successive indices for each
    % unique value fed into it. Returns existing index for a value already seen.
    % Starts by containing a 0 -> 0 map
    
    properties
        map % nT x 2 matrix of [localIndex, globalIndex] pairs
    end
    
    methods
        function this = IndexAccumulator
            this.map = [0, 0];
        end
        
        function globalIndex = index(this, localIndex)
            idxMask = ismember(this.map(:,1), localIndex);
            if any(idxMask) % localIndex already seen
                globalIndex = this.map(idxMask, 2);
            else % new localIndex
                globalIndex = size(this.map, 1);
                this.map = [this.map; localIndex, globalIndex];
            end
        end
    end
end