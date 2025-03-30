classdef incremental_qr < handle
% Implementation of the incremental QR decomposition of a tall matrix
% with an increasing number of columns

    properties
        H   % Matrix to store evolving QR decomposition
        tau % Vector to store scalars for Householder reflections
        k   % Number of columns in current matrix
        n   % Height of matrix
    end
    
    methods
        function obj = incremental_qr(Y,varargin)
            % Initialize the QR decomposition with matrix Y
            [obj.H, obj.tau] = hhqr(Y);
            [obj.n,obj.k] = size(Y);
            if ~isempty(varargin) % Optionally: Set size
                obj.H = [obj.H zeros(obj.n,varargin{1}-obj.k)];
            end
        end
        
        function obj = addcols(obj, Ynew)
            % Add new columns to the existing QR decomposition
            l = size(Ynew,2);

            % Double and copy if necessary
            if obj.k + l > size(obj.H,2)
                obj.H = [obj.H zeros(size(obj.H))];
                obj.tau = [obj.tau;zeros(size(obj.tau))];
                obj.addcols(Ynew);
                return
            end

            Ynew = obj.applyQtfull(Ynew);   % Apply Qt to Ynew
            Ybot = Ynew(obj.k+1:end,:);     % Extract bottom of Ynew

            % Form QR decomposition of bottom
            [Ynew(obj.k+1:end,:),obj.tau(obj.k+1:obj.k+l)]...
                = hhqr(Ybot);

            obj.H(:,obj.k+1:obj.k+l) = Ynew; % Write to buffer
            obj.k = obj.k + l; % Increase size of k
        end
        
        function Qx = applyQ(obj, x)
            % Apply the Q matrix to a vector x
            Qx = matlab.internal.decomposition.applyHouseholder(...
                obj.H, obj.tau, [x;zeros(obj.n-obj.k,...
                size(x,2))], false, obj.k);
        end
        
        function Qtx = applyQt(obj, x)
            % Apply the transpose of the Q matrix to a vector x
            Qtx = apply_Qt(obj.H.obj.tau,x);
            Qtx = Qtx(1:obj.k,:);
        end
        
        function y = projectOut(obj, x)
            % Compute y = (I-Q*Q')*x
            y = apply_Qt(obj.H.obj.tau,x);
            y(1:obj.k,:) = 0;
            y = matlab.internal.decomposition.applyHouseholder(...
                obj.H, obj.tau, y, false, obj.k);
        end
        
        function Q = getQ(obj)
            % Return the Q matrix from the current QR decomposition
            Q = get_Q(obj.H,obj.tau);
        end
        
        function R = getR(obj)
            % Return the R matrix from the current QR decomposition
            R = triu(obj.H);  % Extract the upper triangular part
            R = R(1:obj.k,1:obj.k);
        end
    end
end