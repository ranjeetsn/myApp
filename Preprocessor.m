classdef Preprocessor < handle
    % This class handles the preprocessing of the raw data
    
    properties (Access = private)
        d_PreprocessorType
    end
    
    methods
        
        function obj                = Preprocessor(pp_type)
            obj.SetPreprocessor(pp_type);
        end
        
        function obj                = SetPreprocessor(obj, pp_type)
            obj.d_PreprocessorType =  PreprocessorType.newType(pp_type);
        end
        
        function Preprocessed_data  = PreprocessData(obj, data)
            Preprocessed_data       = obj.d_PreprocessorType.PreprocessData(data);
        end
        
    end
end

