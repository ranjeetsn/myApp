classdef VCV_DiagonalMapper < VarianceCovarianceMapperType
    % This is a concerte class implementing a Variance Covariance Mapper 

    properties
    end
    
    methods
        
        % Maps the variance vector to the diagonal matrix of covariance
        
        function covariance_matrix = Map(obj, variance_vector)
            
          disp(' In the  VCV_DiagonalMapper - Map routine')
          
          covariance_matrix = diag(variance_vector);

        end
        
    end
    
end

