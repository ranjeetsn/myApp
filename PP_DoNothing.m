classdef PP_DoNothing < PreprocessorType
    % This is a concerte class implementing a preprocessor 
    % This preprocessor just returns the input w/o any changes
    properties
    end
    
    methods
        
        function Preprocessed_data = PreprocessData(obj, data)
            
          disp(' In the DoNothing Preprocessor routine')
          Preprocessed_data = data;
        end
        
    end
    
end

