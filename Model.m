classdef Model < handle
    % This class sets the model (static/dynamic)
    
    properties (Access = private)
        d_ModelType
        d_data
        %d_NoOfVariables
%        d_VarianceCovarianceMapper  % Not used anywhere, can be removed
    end
    
    methods

%% Functions in this block initialize the Model (STATIC/DYNAMIC)
        % Initialize the Model object using the model type and the lag
        % For static type the lag should be zero
        function obj                        = Model(md_type)
            obj.m_SetModel(md_type);
        end
        
        function obj                        = m_SetModel(obj, md_type)
            obj.d_ModelType                 = ModelType.newType(md_type);
        end
        
        function obj                        = m_setData(obj, data)
            obj.d_data                      = data;
            obj.d_ModelType.m_setTotalNoOfVariables( size(data.m_getData(),2));
        end
        
        function no_of_variables            = m_getNoOfVariables(obj)
            no_of_variables                 = obj.d_ModelType.m_getTotalNoOfVariables();
        end
        
 %        function obj                       = m_setMapper(obj, vcv_mapper)
 %           obj.d_VarianceCovarianceMapper  = vcv_mapper;
 %        end
        
%% Functions in this block set/configure the DYNAMIC MODEL

        % Function to set the lag and the order of the dynamic model
        function obj                        = m_setLagAndOrderOfModel(obj, lag, order)
            obj.m_setLag(lag);
            obj.m_setOrder(order);
        end

        % Function to set the lag in the model
        function obj                        = m_setLag(obj, lag)
            obj.d_ModelType.m_setLag(lag);
%            obj.m_StackData(obj, lag);
        end
        
        function obj                        = m_setOrder(obj, order)
            obj.d_ModelType.m_setOrder(order);
        end

        % Function to get the lag in the model
        function lag                        = m_getLag(obj)
            lag                             = obj.d_ModelType.m_getLag();
        end
        
        % Function to get the order of the model
        function order                      = m_getOrder(obj)
            order                           = obj.d_ModelType.m_getOrder();
        end
        
        % This function is relevant only if ModelType is DYNAMIC
        function obj                        = m_setNoOfInputsAndOutputs(obj, noOfInputs, noOfOutputs)
            obj.d_ModelType.m_setNoOfInputsAndOutputs(noOfInputs, noOfOutputs);
        end

        % Get the number of Inputs and number of Outputs
        function [noOfInputs, noOfOutputs]  = m_getNoOfInputsAndOutputs(obj)
            [noOfInputs, noOfOutputs]       = obj.d_ModelType.m_getNoOfInputsAndOutputs();
        end
        
        function formatted_data             = m_getFormattedData(obj, data, type)
           formatted_data                   = obj.d_ModelType.m_getFormattedData(data, type); 
        end
        
%%
        % Function to stack the data given in the data obj with the lag in the model
        function lag_stacked_data           = m_StackDataByLag(obj, data)
            lag_stacked_data                = obj.d_ModelType.m_StackDataByLag(data);
        end
        
        % Function to stack the data given in the data obj by the order in the model
        function order_stacked_data         = m_StackDataByOrder(obj, data)
            order_stacked_data              = obj.d_ModelType.m_StackDataByOrder(data);
        end
        
        function [d_min, d_max] = m_getMinAndMaxConstraints(obj)
            [d_min, d_max]      = obj.d_ModelType.m_getMinAndMaxConstraints();
        end
        
        function order = m_ComputeOrderGivenNumOfConstraints(obj, d_val)
            order = obj.d_ModelType.m_ComputeOrderGivenNumOfConstraints(d_val);
        end
        
        % Function to compute the  model
        function model                      = m_computeModel(obj, data_obj, variance_vec, vcv_mapper)
            
           % Here the check should be done on the Noise Model (Noisy/
           % Noisless) and also on the number of arguments
           
            if (nargin == 4)
                model                           = obj.d_ModelType.m_computeModel(data_obj, variance_vec, vcv_mapper);
            elseif (nargin == 2)
                model                           = obj.d_ModelType.m_computeModel(data_obj);
            end
        end
        
        %Function to compute state space parameters
        function m_computeStateSpaceParameters(obj)
            obj.d_ModelType.m_computeStateSpaceParameters();
        end
    end
end
