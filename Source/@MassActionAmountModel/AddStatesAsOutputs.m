function AddStatesAsOutputs(this)
%AddStatesAsOutputs A quick and dirty helper script that adds one output
%   for each state currently in the model. Note: this
%   function is order-dependent, meaning it only matches species already in the
%   model.

full_names = vec(strcat({this.States(1:this.nx).Compartment}, '.', {this.States(1:this.nx).Name}));

for i = 1:numel(full_names)
    this.AddOutput(full_names{i}, full_names{i});
end