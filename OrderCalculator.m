classdef OrderCalculator < handle
    % This class handles the Order Calculations 
    % Internally, it does so by computing the number of constraints
    % using the eigen values. 
    
    
    properties (Access = private)
        
        d_OrderCalculatorType
        d_Order
        d_EigenValues
    end
    
    methods 
        
        function obj                    = OrderCalculator(order_calculator_type)
            obj.m_setOrderCalculator(order_calculator_type);
        end
        
        function obj                    = m_setOrderCalculator(obj, order_calculator_type)
            obj.d_OrderCalculatorType   = OrderCalculatorType.newType(order_calculator_type);
        end

        function obj                    = m_setData(obj, data)
            obj.d_OrderCalculatorType.m_setData(data);
        end
        
        function obj                    = m_setModel(obj, model)
            obj.d_OrderCalculatorType.m_setModel(model);
        end
        
        function obj                    = m_setVarianceModel(obj, variance_model)
            obj.d_OrderCalculatorType.m_setVarianceModel(variance_model);
        end
   
        function obj                    = m_configureOrderCalculator(obj, data, model, variance_model)
            obj.m_setData(data);
            obj.m_setModel(model);
            obj.m_setVarianceModel(variance_model);
        end
                
        function obj                    = m_setOrder(obj, order)
            obj.d_Order                 = order;
        end
        
        function [order, eigen_values]          = m_computeOrder(obj)
            [obj.d_Order,obj.d_EigenValues]     = obj.d_OrderCalculatorType.ComputeOrder();
            order                               = obj.d_Order;
            eigen_values                        = obj.d_EigenValues;
        end
        
        function order  = m_getOrder(obj)
                 order  = obj.d_Order;
        end
        
    end
    
end