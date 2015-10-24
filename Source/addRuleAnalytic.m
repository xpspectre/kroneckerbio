function m = addRuleAnalytic(m, name, target, expression, type)
% Add a rule to a Model.Analytic. Separate rules field simplifies conversion
% between kroneckerbio and SBML formats.
%
% Inputs:
%   m [ Model.Analytic struct ]
%       The model to which the rule will be added
%   name [ string ]
%       Name of rule
%   target [ string ]
%       LHS of rule
%   expression [ string ]
%       RHS of rule
%   type [ {'repeated assignment'} | 'initial assignment' | 'rate' | 'single sub' ]
%       Rule type, taken from SBML spec. 'single sub' is a dummy rule used in
%       NLME fitting.
%
% Outputs:
%   m [ Model.Analytic struct ]
%       The model with the new rule added

% Clean up inputs
if nargin < 5
    type = 'repeated assignment';
end

% Check valid rule type
validRuleTypes = {'repeated assignment', 'initial assignment', 'rate', 'single sub'};
assert(ismember(type, validRuleTypes), 'KroneckerBio:addRuleAnalytic:InvalidRuleType', 'Rule type %s is not a supported rule type', type)

% Increment counter
nz = m.nz + 1;
m.nz = nz;
m.Rules = growRules(m.Rules, m.nz);

% Add item
m.Rules(nz).Name       = name;
m.Rules(nz).Target     = target;
m.Rules(nz).Expression = fixRuleExpression(expression);
m.Rules(nz).Type       = type;

m.Ready = false;
