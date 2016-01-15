function Update(this,k)
% Update model and private components with new parameter values.
%
% Inputs:
%   k [ this.nk array of doubles ]
%       New parameter values

% Clean-up
k = vec(k);
assert(length(k) == this.nk, 'KroneckerBio:Model:Update:InvalidParameterLength', '%g parameters supplied but %g are required.', length(k), this.nk)

% Update private component closure
m = this.m;
m = m.Update(k);
this.m = m;

% Distribute values
if this.nk >= 1
    k = num2cell(k);
    [this.Parameters.Value] = k{:};
end