function [obj, n_obj] = fixObjective(obj, n_con)
% Standardize objective function structures

assert(size(obj,1) == n_con, 'KroneckerBio:objSize', 'obj must have a number of columns equal to numel(con)')

n_obj = size(obj,2);
