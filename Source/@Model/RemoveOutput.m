function RemoveOutput(this, name)
% RemoveOutput Remove an output by name

% Find added instances of this name
ind = strcmp(name, {this.Outputs.Name});
assert(sum(ind)==1, 'KroneckerBio:Model:RemoveOutput:OutputNotFound', 'Output with name %s not found in model', name)

% Remove all mention of this output
this.Outputs(ind,:) = [];
this.ny = this.ny - 1;

this.Ready = false;
