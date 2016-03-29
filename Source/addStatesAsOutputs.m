function m = addStatesAsOutputs(m, include_compartment)
%addStatesAsOutputs A quick and dirty helper script that adds one output
%   for each state currently in the model. Note: this
%   function is order-dependent, meaning it only matches species already in the
%   model.
%
%   m = addStatesAsOutputs(m)
%
% Inputs:
%   m [ model ]
%       The model whose existing states are added as outputs
%   include_compartment [ {true} | false ]
%       Whether to use prepend the compartment to the species names. If false,
%       errors will be thrown on model finalization if multiple outputs refer to
%       species with the same unqualified name.
%
% Outputs:
%   m [ model ]
%       Model with states added as outputs

if nargin < 2
    include_compartment = true;
end

if is(m, 'Model.MassActionAmount')
    for i = 1:m.nx
        if include_compartment
            name = [m.States(i).Compartment '.' m.States(i).Name];
        else
            name = m.States(i).Name;
        end
        
        m = AddOutput(m, name, name);
    end
elseif is(m, 'Model.Analytic')
    for i = 1:m.nx
        if include_compartment
            name = [m.States(i).Compartment '.' m.States(i).Name];
            expression = [quoteIfInvalid(m.States(i).Compartment) '.' quoteIfInvalid(m.States(i).Name)];
        else
            name = m.States(i).Name;
            expression = quoteIfInvalid(m.States(i).Name);
        end
        
        m = AddOutput(m, name, expression);
    end
else
    error('KroneckerBio:AddOutput:m', 'm must be a model')
end

end

function name = quoteIfInvalid(name)
if ~isValidIdentifier(name)
    name = ['"' name '"'];
end
end