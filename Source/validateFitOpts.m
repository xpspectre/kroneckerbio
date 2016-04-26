function isValid = validateFitOpts(opts, models, conditions, objectives)
% Check whether fitting scheme in opts is consistent with the vectors of models,
%   conditions, and objectives passed into it. Used by FitObjective.