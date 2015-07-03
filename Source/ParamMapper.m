classdef (Abstract) ParamMapper < handle
    % Abstract superclass that specifies API of local/condition <-> global/fit
    % parameters mappers
    
    properties (Abstract)
        nT % scalar integer
        paramsMap % nCon x 4 cell matrix of nXi x 1 double vectors
    end
    
    methods (Abstract)
        addCondition(this, condition)
        isSharedParamSpec(this, UseParams) % for addModelOnUniqueParam, useParams [ nk x 1 double vector ]
    end
    
    methods
        function Tlocal = T2Tlocal(this, T)
            % Map overall T vector to local T's in conditions
            % Inputs:
            %    T [ nT x 1 double vector ]
            %       Overall vector of paremeters for fit
            % Outputs:
            %    Tlocal [ nCon x 4 cell matrix of nXi x 1 double vectors ]
            %       Condition-local parameter arrangement
            nCon = size(this.paramsMap, 1);
            Tlocal = cell(nCon,4);
            for i = 1:nCon
                for j = 1:4
                    Tmap_ij = this.paramsMap{i,j};
                    nX = length(Tmap_ij);
                    Tlocal_ij = zeros(size(Tmap_ij));
                    for k = 1:nX
                        Tidx = Tmap_ij(k);
                        if Tidx ~= 0
                            Tlocal_ij(k) = T(Tidx);
                        end
                    end
                    Tlocal{i,j} = Tlocal_ij;
                end
            end
        end
        
        function T = Tlocal2T(this, Tlocal, mode)
            % Map local T's in conditions to overall T vector. Includes optional
            % argument to specify overwriting values in T (multiple conditions
            % sharing a value) or summing (multiple conditions sum gradients)
            % Inputs:
            %    Tlocal [ nCon x 4 cell matrix of nXi x 1 double vectors ]
            %       Condition-local parameter arrangement
            %    mode [ string { 'overwrite' } ]
            %       Controls combination of condition-local parameters allowing:
            %           'overwrite': replace values in baseCell with nonzero values in topCell
            %           'sum': sum values in double vectors
            % Outputs:
            %    T [ nT x 1 double vector ]
            %       Overall vector of paremeters for fit
            if nargin < 3
                mode = 'overwrite';
            end
            
            nCon = size(this.paramsMap, 1);
            T = zeros(this.nT,1);
            for i = 1:nCon
                for j = 1:4
                    Tmap_ij = this.paramsMap{i,j};
                    nX = length(Tmap_ij);
                    Tlocal_ij = Tlocal{i,j};
                    for k = 1:nX
                        Tidx = Tmap_ij(k);
                        if Tidx ~= 0
                            if strcmpi(mode, 'overwrite')
                                T(Tidx) = Tlocal_ij(k);
                            elseif strcmpi(mode, 'sum')
                                T(Tidx) = T(Tidx) + Tlocal_ij(k);
                            end
                        end
                    end
                end
            end
        end
        
    end
    
end