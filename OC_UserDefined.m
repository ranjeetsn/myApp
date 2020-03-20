classdef OC_UserDefined < OrderCalculatorType
    % This is a concerte class implementing a OrderCalculatorType
    % This Class just takes the order information from the user
    
    properties (Access = private)
        d_Order
    end
    
    methods
        
        function obj    = m_setLag(obj, lag)
            obj.d_Lag   = lag;
        end
        
        function obj    = m_setVarianceModel(obj, variance_model)
            obj.d_VarianceModel = variance_model;
        end
        
        function obj    = m_setStackedData(obj, stacked_data)
            obj.d_Data  = stacked_data;
        end
        
        function obj    = m_setOrder(obj, order)
            obj.d_Order = order;
        end
        
        function [order, eigen_values, ret_code] = ComputeOrder(obj)
          disp('In the OC_UserDefined: ComputeOrder routine');
          order         = obj.d_Order;
          disp('Compute and return Eigen Values');
          eigen_values  = [];
          d_val         = obj.d_Lag - obj.d_Order + 1;
          [eigen_values, ret_code]  = SysidUtils.m_ComputeEigenvaluesForGivenD(obj.d_Data, d_val, obj.d_VarianceModel);
        end
        
       
    end
    
end

