classdef PP_SkipAlternate < PreprocessorType
    % This is a concerte class implementing a preprocessor 
    % This preprocessor skips every second element and return the new data
    % object
    properties
    end
    
    methods
        
        function Preprocessed_data = PreprocessData(obj, data)
            
          disp(' In the SkipAlternate Preprocessor routine')
          size(data)
          Preprocessed_data = data(1:10,:);
        end
        
    end
    
end

