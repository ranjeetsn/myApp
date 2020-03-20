classdef ModelType < handle
    % This is an interface/ abstract class for ModelType classes
    % Any class that is implementing a Model needs to inherit this
    % class and define the computeModel method 
    
    properties (Access = protected)
        d_TotalNoOfVariables
        d_Order
    end
    
    
    methods (Abstract)
        
        processed_data = ProcessData(obj, data);
        
        formatted_data  = m_getFormattedData(obj, data, type);
        
        model           =  m_computeModel(obj, data, variance_vec, vcv_mapper);
        [d_min, d_max]  =  m_getMinAndMaxConstraints(obj, l);
        order           =  m_ComputeOrderGivenNumOfConstraints(obj, d_val);
    end
    
    methods 
        
        function obj = m_setOrder(obj, order)
            obj.d_Order = order;
        end
        
        function order = m_getOrder(obj)
            order      = obj.d_Order;
        end
        
        obj = m_setNoOfInputsAndOutputs(obj, noOfInputs, noOfOutputs);
        [noOfInputs, noOfOutputs] = m_getNoOfInputsAndOutputs(obj);
        
       function m_setTotalNoOfVariables(obj, no_of_vars)
            obj.d_TotalNoOfVariables = no_of_vars;
       end
       
       function no_of_vars = m_getTotalNoOfVariables(obj)
             no_of_vars = obj.d_TotalNoOfVariables;
        end
        
    end
    
    methods (Static)
        
        function type  = newType(md_type)
            
            switch md_type
                
                case 'STATIC_MODEL'
                    type = MD_StaticModel;
                    
                case 'DYNAMIC_MODEL'
                    type = MD_DynamicModel;

                otherwise
                    error('Not a valid Model Type');
                    
            end
            
        end
        
    end
    
    
end

