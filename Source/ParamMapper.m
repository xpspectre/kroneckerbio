classdef (Abstract) ParamMapper < handle
    % Abstract superclass that specifies API of local/condition <-> global/fit
    % parameters mappers
    % Order of parameter types: k, s, q, h
    
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
            %    T [ nT x 1 double vector | nT x nT double matrix ]
            %       Overall vector of parameters/gradient or Hessian for fit
            % Outputs:
            %    Tlocal [ nCon x 1 cell vector of 4 x 1 cell vectors of nXi x 1
            %           double vectors | nCon x 1 cell vector of 4 x 4 cell matrices of nXi
            %           x nXj double matrices ]
            %       Condition-local parameters/gradient or Hessian arrangement
            
            % Parameters/gradient (order 1) or Hessian (order 2)
            nTi = size(T,1);
            nTj = size(T,2);
            if nTj == 1
                order = 1;
            elseif nTj > 1 && nTj == nTi
                order = 2;
            else
                error('ParamMapper:T2Tlocal: invalid T dimensions')
            end
            
            nCon = size(this.paramsMap, 1);
            Tlocal = cell(nCon,1);
            
            % Iterate through conditions
            for iCon = 1:nCon
                
                % 4 x 1 for parameters/gradient or 4 x 4 for Hessian
                nXsi = 4;
                if order == 1
                    nXsj = 1;
                else
                    nXsj = 4;
                end
                
                Tmap = this.paramsMap(iCon,:); % condition-specific param map [ 1 x 4 cell vector ]
                Tlocal_iCon = cell(nXsi,nXsj);
                
                % Iterate through parameter types
                for Xsi = 1:nXsi
                    for Xsj = 1:nXsj
                        
                        % nXi x 1 for parameters/gradient or nXi x nXj for Hessian
                        Tmap_i = Tmap{Xsi};
                        nXi = size(Tmap_i,1);
                        if order == 1
                            Tmap_j = 1;
                            nXj = 1;
                        else
                            Tmap_j = Tmap{Xsj};
                            nXj = size(Tmap_j,1);
                        end
                        
                        Tlocal_iCon_ij = zeros(nXi,nXj);
                        
                        % Iterate through positions in Tmap
                        for iX = 1:nXi
                            Tidx_i = Tmap_i(iX);
                            for jX = 1:nXj
                                Tidx_j = Tmap_j(jX);
                                
                                if Tidx_i ~= 0 && Tidx_j ~= 0
                                    Tlocal_iCon_ij(iX,jX) = T(Tidx_i,Tidx_j);
                                end
                            end
                            
                        end
                        
                        Tlocal_iCon{Xsi,Xsj} = Tlocal_iCon_ij;
                        
                    end
                end
                
                Tlocal{iCon,1} = Tlocal_iCon;
                
            end
        end
        
        function T = Tlocal2T(this, Tlocal, mode)
            % Map local T's in conditions to overall T vector. Includes optional
            % argument to specify overwriting values in T (multiple conditions
            % sharing a value) or summing (multiple conditions sum gradients)
            % Inputs:
            %    Tlocal [ nCon x 1 cell vector of 4 x 1 cell vectors of nXi x 1
            %           double vectors | nCon x 1 cell vector of 4 x 4 cell matrices of nXi
            %           x nXj double matrices ]
            %       Condition-local parameters/gradient or Hessian arrangement
            %    mode [ string { 'overwrite' } ]
            %       Controls combination of condition-local parameters allowing:
            %           'overwrite': replace values in T with values extracted
            %               from Tlocal
            %           'sum': sum values in double vectors
            % Outputs:
            %    T [ nT x 1 double vector | nT x nT double matrix ]
            %       Overall vector of parameters/gradient or Hessian for fit
            if nargin < 3
                mode = 'overwrite';
            end
            
            % Parameters/gradient (order 1) or Hessian (order 2)
            nTli = size(Tlocal{1}{1,1},1);
            nTlj = size(Tlocal{1}{1,1},2);
            if nTlj == 1
                order = 1;
            elseif nTlj > 1 && nTlj == nTli
                order = 2;
            else
                error('ParamMapper:Tlocal2T: invalid Tlocal dimensions')
                % TODO: make this error check more thorough
            end
            
%             nCon = size(this.paramsMap, 1);
            nCon = size(Tlocal,1); % Only process the gradients/Hessians passed to it
            if order == 1
                T = zeros(this.nT,1);
            else
                T = zeros(this.nT,this.nT);
            end
            
            % Iterate through conditions
            for iCon = 1:nCon
                
                % 4 x 1 for parameters/gradient or 4 x 4 for Hessian
                nXsi = 4;
                if order == 1
                    nXsj = 1;
                else
                    nXsj = 4;
                end
                
                Tmap = this.paramsMap(iCon,:); % condition-specific param map [ 1 x 4 cell vector ]
                Tlocal_iCon = Tlocal{iCon};
                
                % Make Tlocal_iCon a column vector for order 1 for safety
                if order == 1
                    Tlocal_iCon = vec(Tlocal_iCon);
                end
                
                % Iterate through parameter types
                for Xsi = 1:nXsi
                    for Xsj = 1:nXsj
                        
                        % nXi x 1 for parameters/gradient or nXi x nXj for Hessian 
                        Tmap_i = Tmap{Xsi};
                        nXi = size(Tmap_i,1);
                        if order == 1
                            Tmap_j = 1;
                            nXj = 1;
                        else
                            Tmap_j = Tmap{Xsj};
                            nXj = size(Tmap_j,1);
                        end
                        
                        Tlocal_iCon_ij = Tlocal_iCon{Xsi,Xsj};
                        
                        % Iterate through positions
                        for iX = 1:nXi
                            Tidx_i = Tmap_i(iX);
                            for jX = 1:nXj
                                Tidx_j = Tmap_j(jX);
                                
                                if Tidx_i ~= 0 && Tidx_j ~= 0
                                    if strcmpi(mode, 'overwrite')
                                        T(Tidx_i,Tidx_j) = Tlocal_iCon_ij(iX,jX);
                                    elseif strcmpi(mode, 'sum')
                                        T(Tidx_i,Tidx_j) = T(Tidx_i,Tidx_j) + Tlocal_iCon_ij(iX,jX);
                                    end
                                end
                            end
                        end
                        
                    end
                end
                
            end
            
        end
        
    end
    
end