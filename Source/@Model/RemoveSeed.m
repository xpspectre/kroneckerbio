function RemoveSeed(this, name)
% RemoveSeed Remove a seed by name

% Find added instances of this name
ind = strcmp(name, {this.Seeds.Name});
assert(sum(ind)==1, 'KroneckerBio:Model:RemoveSeed:SeedNotFound', 'Seed with name %s not found in model', name)

% Remove all mention of this seed
this.Seeds(ind,:) = [];
this.ns = this.ns - 1;

this.Ready = false;