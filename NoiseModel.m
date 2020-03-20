classdef NoiseModel < handle
    %This singleton implementation for the Noise Model provides 
    %a place to store the noise model type and single point of access
    
    properties (Access = private)
        d_NoiseModelType
    end
    
    methods(Access=private)
        
        function obj = NoiseModel()
            obj.d_NoiseModelType = NoiseModelEnums.UNDEFINED;
        end
        
    end
    
    methods (Static)
        
        function obj = getInstance()
            
            persistent d_instance;
            
            if (isempty(d_instance))
                obj = NoiseModel();
                d_instance = obj;
            else
                obj = d_instance;
            end
        end
        
    end
    
    methods
         
        function obj                    = m_setNoiseModel(obj, noise_model_type)
                obj.d_NoiseModelType    = noise_model_type;
        end
        
        function noise_model            = m_getNoiseModel(obj)
            noise_model                 = obj.d_NoiseModelType;
        end
        
    end
end

