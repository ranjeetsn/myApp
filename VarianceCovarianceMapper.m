classdef VarianceCovarianceMapper < handle
    % This class handles the preprocessing of the raw data
    
    properties (Access = private)
        d_VarianceCovarianceMapperType
    end
    
    methods
        
        function obj                = VarianceCovarianceMapper(vcv_mapper_type)
            obj.m_setVarianceCovarianceMapper(vcv_mapper_type);
        end
        
        function type = m_getType(obj)
            type = obj.d_VarianceCovarianceMapperType;
        end
        
        function obj                           =  m_setVarianceCovarianceMapper(obj, vcv_mapper_type)
            obj.d_VarianceCovarianceMapperType =  VarianceCovarianceMapperType.newType(vcv_mapper_type);
        end
        
        % This function is to be called only when
        % VarianceCovarianceMapperType is DipcaMapper
        function obj    = m_setNoOfVariablesAndLag(obj, no_of_variables, lag)
            obj.d_VarianceCovarianceMapperType.m_SetNoOfVariables(no_of_variables);
            obj.d_VarianceCovarianceMapperType.m_SetLag(lag);
        end
        
        function obj    = m_setLag(obj, lag)
            obj.d_VarianceCovarianceMapperType.m_SetLag(lag);
        end
        
        % This function is to be called only when 
        % VarianceCovarianceMapperType is GenericMapper
        function obj = m_setRowIndicesAndColumnIndices(obj, row_ind, col_ind)
            obj.d_VarianceCovarianceMapperType.m_setRowIndexAndColumnIndex(row_ind, col_ind);
        end
         
        function covariance_matrix  = m_mapVarianceToCovariance(obj, variance_vector)
            covariance_matrix       = obj.d_VarianceCovarianceMapperType.Map(variance_vector);
        end
        
    end
end

