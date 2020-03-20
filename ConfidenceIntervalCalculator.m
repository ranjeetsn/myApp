classdef ConfidenceIntervalCalculator < handle

    
    properties
        d_Data
        d_NumberOfBootstraps
        d_SizeOfBootstraps
        d_AlphaValue
    end
    
    methods
        
        function obj                    = m_setData(obj, data)
            obj.d_Data                  = data;
        end
        
        function obj                    = m_setNumberOfBootstraps(obj, no_of_bootstraps)
            obj.d_NumberOfBootstraps    = no_of_bootstraps;
        end
        
        function obj                    = m_setSizeOfBootstraps(obj, size_of_bootstraps)
            % Check if size of bootstrap is greater than the total number
            % of rows and throw an error
            obj.d_SizeOfBootstraps      = size_of_bootstraps;
        end
        
        function obj                    = m_setAlphaValue(obj, alpha_value)
            obj.d_AlphaValue            = alpha_value;
        end
        
        function confidence_intervals   = m_computeConfidenceIntervals(obj, order, n_outputs, lag, variance_vec, vcv_mapper)
            
            model       = [];
            confidence_intervals = [];
            data        = obj.d_Data.m_getData();
            total_cols  = size(data, 2);
            nvar        = size(data,1);
            
            for i = 1:obj.d_NumberOfBootstraps
                
                col_indices     = randperm(total_cols, obj.d_SizeOfBootstraps);
                selected_rows   = data(:,col_indices);
                
                if (nargin == 6)
                    temp            = SysidUtils.m_computeModelForStackedData(selected_rows, order, ...
                                                                              n_outputs, lag, ...
                                                                              variance_vec,  vcv_mapper);
                elseif (nargin == 4)
                     temp            = SysidUtils.m_computeModelForStackedData(selected_rows, order, ...
                                                                              n_outputs, lag);
                                                                              
                end
                
                temp            = temp/temp(nvar,1);
                model           = [model temp];
                
            end
            mean_model = mean(model, 2)
            std_model  = std(model, 0, 2);
            
            for i = 1:nvar-1
                interval = norminv([obj.d_AlphaValue/200 1-obj.d_AlphaValue/200], mean_model(i), std_model(i));
                 confidence_intervals = [ confidence_intervals; interval];
            end
            confidence_intervals = [ confidence_intervals; 1 1];
           
        end
        
    end
end

