function RemoveCompartment(this, name)
% RemoveCompartment Remove a compartment by name

% Find added instances of this name
ind = strcmp(name, {this.Compartments.Name});
assert(sum(ind)==1, 'KroneckerBio:Model:RemoveCompartment:CompartmentNotFound', 'Compartment with name %s not found in model', name)

% Remove all mention of this compartment
this.Compartments(ind,:) = [];
this.nv = this.nv - 1;

this.Ready = false;