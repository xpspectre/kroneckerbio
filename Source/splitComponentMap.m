function componentMaps = splitComponentMap(componentMap, n)
% Split component map by condition into n ~equal chunks. componentMap has a row
%   for each obj fun, so combine all obj funs for a particular condition
%   together.
conditions = unique(componentMap(:,2));
conditionChunks = splitList(conditions, n);
componentMaps = cell(n,1);
for in = 1:n
    conditions = conditionChunks{in};
    components = zeros(0,3);
    for iCon = 1:length(conditions)
        condition = conditions(iCon);
        components = [components; componentMap(componentMap(:,2)==condition,:)]; % resizing shouldn't be a problem for something this small
    end
    componentMaps{in} = components;
end
end