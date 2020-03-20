classdef VR_KnownVariance < VarianceModelType
    % This is a concerte class implementing a Model 
    % This Model just returns the input w/o any changes
    properties
    end
    
    methods
        
        function covariance = m_ComputeCovariance(obj, data, num_of_constraints)
            
          disp(' In the VR_KnownVariance : m_ComputeCovariance')
          
          covariance = [];
          
          % set the variances to the known variances taken from user in
          % data
        end
        
    end
    
end

