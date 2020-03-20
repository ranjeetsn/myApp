classdef VarianceCovarianceMapperType < handle
    % This is an interface/ abstract class for mapping classes
    % Any class that is implementing a mapping from variance to covariance
    %needs to inherit this
    
    properties
    end
    
    
    methods (Abstract)
        covariance_matrix = Map(obj, variance_vector);
    end
    
    
    methods (Static)
        
        function type = newType(vcv_mapper_type)
            
            switch vcv_mapper_type
                
                case 'DIAGONAL_MAPPER'
                    type = VCV_DiagonalMapper;
                    
                case 'DIPCA_MAPPER'
                    type = VCV_DipcaMapper;
                    
                case 'GENERIC_MAPPER'
                    type = VCV_GenericMapper;
                    
                otherwise
                    error('Not a valid Preprocessor Type');   
            end
            
        end
        
    end

end

