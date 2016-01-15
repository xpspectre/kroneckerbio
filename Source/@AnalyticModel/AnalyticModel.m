classdef AnalyticModel < Model
    
    methods (Access = protected)
        qualifyNames(this)
        calculateDerivatives(this, opts)
    end
    
    methods
        function this = AnalyticModel(name)
            if nargin < 1
                name = [];
            end
            this@Model(name);
        end
    end
end