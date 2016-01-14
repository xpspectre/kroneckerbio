function AddState(this, name, compartment, ic)
%AddState Add a state species to a model. The state's initial
%condition can be a scalar value or an arbitrary expression of seeds.
%
%   Inputs
%   m: [ model struct scalar ]
%       The model to which the state will be added
%   name: [ string ]
%       A name for the state. This is the name by which reactions will
%       refer to it.
%   compartment: [ string ]
%       The name of the compartment to which it will be added
%   ic: [ string | nonnegative scalar {0} ]
%       String expression for the initial condition in terms of seeds or
%       nonnegative scalar initial amount.

import Validate.*

% Clean up inputs
if nargin < 4
    ic = [];
end

% Increment counter
nx = this.nx + 1;
this.nx = nx;
this.growStates;

% Add item
this.States(nx).Name         = Validate.SpeciesName(name);
this.States(nx).Compartment  = Validate.CompartmentName(compartment);
this.States(nx).InitialValue = Validate.StateInitialCondition(ic);

this.Ready = false;
