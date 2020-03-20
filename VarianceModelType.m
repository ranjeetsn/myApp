classdef VarianceModelType
    % This is an interface/ abstract class for covariance computation classes
    % Any class that is implementing a Covariance calculation needs to inherit this
    % class and define the ComputeCovariance method with Dataobj object as
    % input and covariance matrix as output
    
    properties
    end
    
    
    methods (Abstract)
        covariance = m_ComputeCovariance(obj, data, no_of_constraints);
    end
    
    
    methods (Static)
        
        function type = newType(vr_type)
            
            switch vr_type
                
                case 'KNOWN'
                    type = VR_KnownVariance;
                    
                case 'UNKNOWN'
                    type = VR_UnknownVariance;
                    
                otherwise
                    error('Not a valid Model Type');
                    
            end
            
        end
        
    end
    
    
end

