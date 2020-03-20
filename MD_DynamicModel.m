classdef MD_DynamicModel < ModelType
    
    % This is a concerte class implementing a Model 
    % The Dynamic Model implements a stacking routine 
    % for converting a dynamic model to a static model
    % Also, to compute the model given data and error covariance
    % information
    
    properties (Access = private)
        d_Lag           = 0 
%        d_Order         = 0
        d_noOfInputs    = 0
        d_noOfOutputs   = 0
        d_constraintMatrix 
    end
    
    methods
%% Methods to set the properties of the class        

        function obj        = m_setLag(obj, lag)
            obj.d_Lag       = lag;
        end
        
%        function obj        = m_setOrder(obj, order)
%            obj.d_Order     = order;
%        end
        
        function lag        = m_getLag(obj)
            lag             = obj.d_Lag;
        end
        
%         function order      = m_getOrder(obj)
%             order           = obj.d_Order;
%         end
        
        function obj = m_setNoOfInputsAndOutputs(obj, noOfInputs, noOfOutputs)
        % Include a sanity check to see if sum of noOfInputs and noOfOutputs add up to dimensions of data matrix
            obj.d_noOfInputs    = noOfInputs;
            obj.d_noOfOutputs   = noOfOutputs;
            obj.d_TotalNoOfVariables = noOfInputs + noOfOutputs;
        end

        % Get the number of Inputs and number of Outputs
        function [noOfInputs, noOfOutputs] = m_getNoOfInputsAndOutputs(obj)
            noOfInputs                     = obj.d_noOfInputs;
            noOfOutputs                    = obj.d_noOfOutputs;
        end
        
        
%% Methods to process data        

        function formatted_data  = m_getFormattedData(obj, data, type)
            if (nargin == 2)
                formatted_data            = SysidUtils.GenerateStackedMatrix(data.m_getData(), obj.d_Lag);
            elseif (nargin == 3)
                if (type == StackingTypeEnums.STACK_BY_ORDER)
                    formatted_data            = SysidUtils.GenerateStackedMatrix(data.m_getData(), obj.d_Order);
                else 
                    formatted_data            = SysidUtils.GenerateStackedMatrix(data.m_getData(), obj.d_Lag);
                end
            end
            
        end
        
        function stacked_data       = ProcessData(obj, data)
            stacked_data            = m_StackDataByLag(obj, data);
        end
        
        function stacked_data       = m_StackDataByLag(obj, data)
            disp(' In the MD_DynamicModel: m_StackDataByLag routine')
            stacked_data            = SysidUtils.GenerateStackedMatrix(data.m_getData(), obj.d_Lag);
            
        end
        
        function stacked_data       = m_StackDataByOrder(obj, data)
            disp(' In the MD_DynamicModel: m_StackDataByOrder routine')
            stacked_data            = SysidUtils.GenerateStackedMatrix(data.m_getData(), obj.d_Order);
        end
        
        function [d_min, d_max]     = m_getMinAndMaxConstraints(obj)
            d_min                   = 2;
            d_max                   = obj.d_noOfOutputs * (obj.d_Lag + 1) - 1;
        end
        
        function order           =  m_ComputeOrderGivenNumOfConstraints(obj, d_val)
            order = obj.d_noOfOutputs * ( obj.d_Lag + 1) - d_val;
        end

%% Method that computes the model given the error variances       

        function model              = m_computeModel(obj, input_data, variance_vec, vcv_mapper)
            if (nargin == 4)
                model                   = obj.m_computeModelStackedByOrder(input_data, variance_vec, vcv_mapper);
            else
                model                   = obj.m_computeModelStackedByOrder(input_data);
            end
        end

        function model = m_computeModelStackedByLag(obj, input_data, variance_vec, vcv_mapper)
            
            noise_model             = NoiseModel.getInstance();
            stacked_data            = obj.m_StackDataByLag(input_data);
            
            if (noise_model.m_getNoiseModel() == NoiseModelEnums.NOISY)    
               model                   = SysidUtils.m_computeModelForStackedData(stacked_data, obj.d_Order, ...
                                                                                  obj.d_noOfOutputs, obj.d_Lag, ...
                                                                                  variance_vec, vcv_mapper);
            else
                % Compute constraints for the Noiseless case 
                 model                   = SysidUtils.m_computeModelForStackedData(stacked_data, obj.d_Order, ...
                                                                                  obj.d_noOfOutputs, obj.d_Lag);
                                                                                  
            end
            obj.d_constraintMatrix  = model';
        end
        

        function model              = m_computeModelStackedByOrder(obj, input_data, variance_vec, vcv_mapper)
           
            noise_model             = NoiseModel.getInstance();
            stacked_data            = obj.m_StackDataByOrder(input_data);
            
            if (noise_model.m_getNoiseModel() == NoiseModelEnums.NOISY)
                    
                model                   = SysidUtils.m_computeModelForStackedData(stacked_data, obj.d_Order, ...
                                                                                  obj.d_noOfOutputs, obj.d_Order, ...
                                                                                  variance_vec, vcv_mapper);
                
            else
                % Compute constraints for the Noiseless case 
                model                   = SysidUtils.m_computeModelForStackedData(stacked_data, obj.d_Order, ...
                                                                                  obj.d_noOfOutputs, obj.d_Order);
                
            end
            obj.d_constraintMatrix  = model';
        end   
        
%% Methods handling the state space parameters
        function obj                = m_computeStateSpaceParameters(obj)
            
            num_of_output_vars      = obj.d_noOfOutputs * (obj.d_Lag + 1);
            gamma_s_perpn_transpose = obj.d_constraintMatrix(:,1:num_of_output_vars);
            negative_gamma_s_perpn_transpose_Hs = obj.d_constraintMatrix(:,(num_of_output_vars + 1):end);
            
            m = obj.d_noOfOutputs;
            s = obj.d_Lag;
            n = obj.d_Order;
            
            gamma_s     = null(gamma_s_perpn_transpose);
            gamma_x     = [];
            gamma_y     = [];
            
            for i = 1:obj.d_Lag
                gamma_x = [gamma_x ; gamma_s( (m * (i-1) + 1):(m * s),: )];
            end
            
            for i = 1:obj.d_Lag
                gamma_y = [gamma_y ; gamma_s( (m * i + 1):(m * (s + 1)),: )];
            end
            
            AT = inv( gamma_y' * gamma_y ) * (gamma_y' * gamma_x);
            
            CT = gamma_s( (m * s + 1):(m * (s + 1)), :);
            
            y_coeff = charpoly(AT)
            
            gamma_s_1_perpn_transpose = [];
            
            for i = 1:(s + 1)
                gamma_s_1_perpn_transpose = [ gamma_s_1_perpn_transpose; y_coeff(i) * eye(m) ];
            end
            
            Hs = -(pinv(gamma_s_perpn_transpose)) * (negative_gamma_s_perpn_transpose_Hs);
            
            u_coeff = gamma_s_1_perpn_transpose'* Hs
            
            PT = obj.d_constraintMatrix;
            
            phi = - PT(:, 1:(m * (s + 1)) );
            
            psi = PT(:,((m * (s + 1)) + 1):end);
            
            stacked_phi = phi;
            
            padding_matrix = zeros ( (m * (s+1) - n), m);
            
            new_stack_row = stacked_phi;
            
            for i = 1:s
                new_mat = [padding_matrix new_stack_row];
                new_stack_row = new_mat( :, 1:(m * (s+1)) );
                stacked_phi = [ new_stack_row; stacked_phi];
            end

            h =  pinv(stacked_phi'*stacked_phi) * (stacked_phi'*psi(:));
            
            gamma_m1_ms = gamma_s( m+1:m*(s + 1), :);
            
            left_mat = [ gamma_m1_ms zeros( m * s, m);
                        zeros(m, n)  eye(m,m)];
           
            b_d_estimate = pinv(left_mat) * h        
            
        end
    end
    
end

