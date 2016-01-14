function RemoveParameter(this, name)
% RemoveParameter Remove a parameter by name

% Find added instances of this name
ind = strcmp(name, {this.Parameters(1:this.nk).Name});
assert(sum(ind)==1, 'KroneckerBio:Model:RemoveParameter:ParameterNotFound', 'Parameter with name %s not found in model', name)

% Remove all mention of this parameter
this.Parameters(ind,:) = [];
this.nk = this.nk - 1;

this.Ready = false;
