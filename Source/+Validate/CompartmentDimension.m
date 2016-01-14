function dimension = CompartmentDimension(dimension)
% Standardize the compartment dimension as 0, 1, 2, or 3

assert(isnumeric(dimension) && isscalar(dimension) && any(dimension == 0:3), ...
    'KroneckerBio:Validator:CompartmentDimension', 'Compartment dimension must be 0, 1, 2, or 3')
