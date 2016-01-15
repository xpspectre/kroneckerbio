classdef MassActionAmountModel < Model
    
    methods (Access = protected)
        growReactions(this) % overrides Model's use of Rate with Parameter in Reactions
        qualifyNames(this)
        calculateDerivatives(this)
    end
    
    methods
        function this = MassActionAmountModel(name)
            if nargin < 1
                name = [];
            end
            this@Model(name);
        end
    end
end