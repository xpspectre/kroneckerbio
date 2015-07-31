classdef ParamMapperMultiModelType < ParamMapper
   % Specifies local/condition <-> global/fit parameters mapping with multiple
   % model types where within each parameter type, indices refer to unique
   % parameters. Needed to line up disparate parameters in different models
   % (e.g. same rate constant in 2 model types with different topologies)
   
   properties
       nT = 0;
       paramsMap = cell(0,4);
       paramsShared = cell(0,4);
   end
   
   methods
       function this = ParamMapperMultiModelType(conditions)
           if nargin < 1
               conditions = [];
           end
           
           for i = 1:length(conditions)
               this.addCondition(conditions(i));
           end
       end
       
       function addCondition(this, condition)
           assert(~isempty(condition.ParamsSpec), 'ParamMapperMultiModelType:addCondition: ParamsSpec not found');
           assert(iscell(condition.ParamsSpec), 'ParamMapperMultiModelType:addCondition: ParamsSpec not a cell array');
           assert(length(condition.ParamsSpec) == 4, 'ParamMapperMultiModelType:addCondition: ParamsSpec has wrong size');
           
           nCon = size(this.paramsShared, 1) + 1;
           for i = 1:4
               this.paramsShared{nCon,i} = condition.ParamsSpec{i};
           end
           
           this.updateParamsMap;
       end
       
       function sharedIdx = isSharedParamSpec(this, UseParams)
           % Tests whether condition's UseParams matches an existing one -
           % used by addModelOnUniqueParam to determine whether a new model
           % needs to be created.
           % Note: Doesn't handle corner case when different model types have
           % the same number of k and try to fit the exact same values to them
           % (don't know why you would ever do this, though)
           sharedIdx = 0;
           nCon = size(this.paramsShared, 1);
           for i = 1:nCon
               if length(UseParams) ~= length(this.paramsShared{i,1}) % skip if different model type and has different number of k
                   continue
               end
               if all(UseParams == this.paramsShared{i,1});
                   sharedIdx = i;
                   return
               end
           end
       end
       
       function updateParamsMap(this)
           % Create paramsMap and calculate nT from basic paramsShared cell
           % array from individual conditions.
           Toffset = 0;
           paramsMap_ = cell(size(this.paramsShared));
           
           for i = 1:4
               XShared = this.paramsShared(:,i);
               XMap = paramsShared2paramsMap(XShared);
               [XMap, Toffset] = addOffsetToNonzeroCellArray(XMap, Toffset);
               paramsMap_(:,i) = XMap;
           end
           
           this.paramsMap = paramsMap_;
           this.nT = Toffset;
       end
       
   end
end

function paramsMap = paramsShared2paramsMap(paramsShared)
% Turn individual condition parameter spec into overall mapping
% Input:
%    paramsShared [ nCon x 1 cell vector of nXi x 1 double vectors ]
% Output:
%    paramsMap [ nCon x 1 cell vector of nXi x 1 double vectors ]

ia = IndexAccumulator;
paramsMap = cell(size(paramsShared));
nCon = length(paramsShared);
for i = 1:nCon
    paramsShared_i = paramsShared{i};
    paramsMap_i = zeros(size(paramsShared_i));
    for j = 1:length(paramsShared_i)
        paramsMap_i(j) = ia.index(paramsShared_i(j));
    end
    paramsMap{i} = paramsMap_i;
end

end

function [cellOut, offsetOut] = addOffsetToNonzeroCellArray(cellIn, offsetIn)
% Add offset to all nonzero values in cell array of double vectors
% Optionally specify input (default 0) and output offset, where offset is the
% max value found in cellOut
% Input:
%   cellIn [ nCon x 1 cell vector of nXi x 1 double vectors ]
%   offsetIn [ scalar integer ]
% Output:
%   cellOut [ nCon x 1 cell vector of nXi x 1 double vectors ]
%   offsetOut [ scalar integer ]
if nargin < 2
    offsetIn = 0;
end

maxIdx = 0;
cellOut = cell(size(cellIn));
nCon = length(cellIn);
for i = 1:nCon;
    if ~isempty(cellIn{i}) && max(cellIn{i}) > maxIdx
        maxIdx = max(cellIn{i});
    end
    cellOut{i} = (cellIn{i} + offsetIn) .* logical(cellIn{i}); % selecting with logical restores 0's filled in by offset
end

offsetOut = offsetIn + maxIdx;
end