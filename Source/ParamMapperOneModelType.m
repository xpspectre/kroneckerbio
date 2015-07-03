classdef ParamMapperOneModelType < ParamMapper
   % Specifies local/condition <-> global/fit parameters mapping with 1 type of
   % model
   
   properties
       nT = 0;
       paramsMap = cell(0,4);
       paramsShared = cell(0,4);
   end
   
   methods
       function this = ParamMapperOneModelType(conditions)
           if nargin < 1
               conditions = [];
           end
           
           for i = 1:length(conditions)
               this.addCondition(conditions(i));
           end
       end
       
       function addCondition(this, condition)
           assert(~isempty(condition.ParamsSpec), 'ParamMapperOneModelType:addCondition: ParamsSpec not found');
           assert(iscell(condition.ParamsSpec), 'ParamMapperOneModelType:addCondition: ParamsSpec not a cell array');
           assert(length(condition.ParamsSpec) == 4, 'ParamMapperOneModelType:addCondition: ParamsSpec has wrong size');
           
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
           sharedIdx = 0;
           nCon = size(this.paramsShared, 1);
           for i = 1:nCon
               if all(UseParams == this.paramsShared{i,1});
                   sharedIdx = i;
               end
           end
       end
       
       function updateParamsMap(this)
           % Create paramsMap and calculate nT from basic paramsShared cell
           % array from individual conditions.
           Toffset = 0;
           nCon = size(this.paramsShared, 1);
           paramsMap_ = cell(size(this.paramsShared));
           
           for i = 1:4
               XShared = [this.paramsShared{:,i}]; % combine and reshape into a matrix since all conditions share the same model param spec
               XMap = paramsShared2paramsMap(XShared);
               [XMap, Toffset] = addOffsetToNonzero(XMap, Toffset);
               nX = size(XMap, 1);
               paramsMap_(:,i) = mat2cell(XMap, nX, ones(1,nCon))'; % split and reshape into expected cell matrix
           end
           
           this.paramsMap = paramsMap_;
           this.nT = Toffset;
       end
       
   end
end

function paramsMap = paramsShared2paramsMap(paramsShared)
% Turn individual condition parameter spec into overall mapping
% Input:
%    paramsShared [ nX x nCon double matrix ]
% Output:
%    paramsMap: [ nX x nCon double matrix ]
[nX, nCon] = size(paramsShared);
paramsMap = zeros(nX, nCon);
Toffset = 0;
for i = 1:nX
    XUse = paramsShared(i,:);
    XNonZero = logical(XUse);
    XUseNonzero = XUse(XNonZero);
    [c,~,ic] = unique(XUseNonzero, 'stable');
    
    Tidx = ic' + Toffset;
    
    paramsMap(i,XNonZero) = Tidx;
    Toffset = Toffset + nnz(c);
end
end

function [matOut, offsetOut] = addOffsetToNonzero(matIn, offsetIn)
% Add offset to all nonzero values in matrix.
% Optionally specify input (default 0) and output offset, where additional
% offset is the max value found in matOut
% Inputs:
%   matIn [ nX x nCon double matrix ]
%   offsetIn [ scalar integer ]
% Outputs:
%   matOut [ nX x nCon double matrix ]
%   offsetOut [ scalar integer ]
if nargin < 2
    offsetIn = 0;
end

matOut = (matIn + offsetIn) .* logical(matIn); % selecting with logical restores 0's filled in by offset

if isempty(matIn)
    offsetOut = offsetIn;
else
    offsetOut = offsetIn + max(max(matIn));
end
end
