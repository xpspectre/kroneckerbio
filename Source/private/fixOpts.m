function newOpts = fixOpts(m, opts)

nm = length(m);
UseParams = opts.UseParams;

newOpts(nm,1) = opts;
for i = 1:nm
    tempOpts = opts;
    
    % Get varying params
    tempOpts.UseParams = logical(UseParams(:,i));
    tempOpts.UseSeeds = opts.UseSeeds(:,i);
    tempOpts.UseInputControls = opts.UseInputControls(i);
    tempOpts.UseDoseControls = opts.UseDoseControls(i);
    
    % Get bounds
    tempOpts.LowerBound = kVec2k(opts.LowerBound, UseParams, i);
    tempOpts.UpperBound = kVec2k(opts.UpperBound, UseParams, i);
    
    % Fix tolerances
    tempOpts.AbsTol = opts.AbsTol(i);
    
    % Misc
    tempOpts.ObjWeights = opts.ObjWeights(i,i);
    tempOpts.continuous = opts.continuous(i);
    tempOpts.complex = opts.complex(i);
    tempOpts.tGet = opts.tGet{i};
    
    newOpts(i) = tempOpts;
end