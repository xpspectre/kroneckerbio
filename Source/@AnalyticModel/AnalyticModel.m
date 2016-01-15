classdef AnalyticModel < Model
    
    methods (Access = protected)
        qualifyNames(this)
        calculateDerivatives(this)
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