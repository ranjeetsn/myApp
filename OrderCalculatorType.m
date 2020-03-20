classdef OrderCalculatorType < handle
    % This is an interface/ abstract class for Order determination classes
    % Any class that is implementing a OrderCalculator needs to inherit this
    % class and define the ComputeOrder method
    
    properties (Access = protected)
        d_Data
        d_Model
        d_VarianceModel
    end
    
    methods (Abstract)
        [order, eig_vals, ret_code] = ComputeOrder(obj);
        m_setOrder(order);
    end
    
    methods (Static)
        
        function type  = newType(order_calculator_type)
            
            switch order_calculator_type
                
                case 'USER_DEFINED'
                    type = OC_UserDefined;
                    
                case 'HYPOTHESIS_TESTING'
                    type = OC_HypothesisTesting;
                    
                case 'EIGEN_PLOTS'
                    type = OC_EigenPlots;
                  
                otherwise
                    error('Not a valid OrderCalculator Type');   
                    
            end 
        end
        
    end
    
    methods
        
        function obj                    = m_setData(obj, data)
            obj.d_Data                  = data;
        end
        
        function obj                    = m_setModel(obj, model)
            obj.d_Model                 = model;
        end
        
        function obj                    = m_setVarianceModel(obj, variance_model)
            obj.d_VarianceModel         = variance_model;
        end
       
    end
end

