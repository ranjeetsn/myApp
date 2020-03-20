classdef MD_StaticModel < ModelType
    % This is a concerte class implementing a Model 
    % This Model just returns the input w/o any changes
    properties
    end
    
    methods
        
        function processed_data = ProcessData(obj, data)
            
          disp(' In the MD_StaticModel: ProcessData routine')
          processed_data        = data;
        end
        
        
        function formatted_data = m_getFormattedData(obj, data, ~)
            formatted_data      = data.m_getData()';
        end
 
        function [d_min, d_max] = m_getMinAndMaxConstraints(obj)
            d_min               = 2;
            d_max               = obj.d_TotalNoOfVariables - 1;
        end
        
        function order          = m_ComputeOrderGivenNumOfConstraints(obj, d_val)
                order           = obj.d_TotalNoOfVariables - d_val;
        end
        
        function model          = m_computeModel(obj, input_data, variance_vec, vcv_mapper)
            
            data                = input_data.m_getData();
            data                = data';
            
            no_of_constraints   = obj.d_TotalNoOfVariables - obj.d_Order;
            noise_model         = NoiseModel.getInstance();
            
            if (noise_model.m_getNoiseModel() == NoiseModelEnums.NOISY && nargin == 4) 
                Qe              = vcv_mapper.m_mapVarianceToCovariance(variance_vec);    
                model           = SysidUtils.m_computeConstraints( data, no_of_constraints, Qe);
            else
                % Compute constraints for the Noiseless case 
                model           = SysidUtils.m_computeConstraints( data, no_of_constraints);
            end
        end
        
        
    end
    
end

