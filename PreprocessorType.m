classdef PreprocessorType
    % This is an interface/ abstract class for preprocessing classes
    % Any class that is implementing a preprocessor needs to inherit this
    % class and define the PreprocessData method with Dataobj object as
    % input and Dataobj object as output
    
    properties
    end
    
    
    methods (Abstract)
        Preprocessed_data = PreprocessData(obj, data);
    end
    
    
    methods (Static)
        
        function type = newType(pp_type)
            
            switch pp_type
                
                case 'NULL'
                    type = PP_DoNothing;
                    
                case 'SKIP_ALTERNATE'
                    type = PP_SkipAlternate;
                    
                    
                otherwise
                    error('Not a valid Preprocessor Type');   
            end
            
        end
        
    end

end

