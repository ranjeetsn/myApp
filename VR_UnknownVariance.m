classdef VR_UnknownVariance < VarianceModelType
    
    % This is a concerte class implementing a Model when error variances
    % are unknown
    
    properties
        d_NoOfConstraints = 0  % Applicable for dynamic models. Need for computing the constraint matrix 
        d_VarianceModel
    end
    
    methods
        
        
        function obj                = m_setNoOfConstraints(obj, no_of_constraints)
            obj.d_NoOfConstraints   = no_of_constraints;
        end
        
        
        function covariance = m_ComputeCovariance(obj)
            obj.d_VarianceModel = VarianceModel(
          disp(' In  VR_UnknownVariance : m_ComputeCovariance')
          
          % should impose the DIPCA constraint here thro mapper and incorporate into
          % IPCA
        %  covariance = [];
          
       %   data = dataobj.d_data;
          %
             
          %Compute the variances using the IPCA algo and set the variances
          %the mapping function should come here
        end     
      
      
    end
end
