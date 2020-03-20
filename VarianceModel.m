classdef VarianceModel < handle
%%  This class handles the variance computation. 
%   Variance Model can be (KNOWN/ UNKNOWN)
%   Error variances can set set if they are KNOWN
%   If error variances are UNKNOWN set up the 
%   Variance - Covariance Mapper & the Variance Solver

    
    properties (Access = private)
        d_VarianceModelType
        d_errorVariances
        d_defaultErrorVariances
        d_VarianceComputationStatus
    end
    
    properties
        d_varianceCovarianceMapper
        d_varianceSolver
    end
%%    
    methods
%%  The functions in this block initialize the basic  variance model with type(KNOWN/UNKNOWN) 
%   and Error Variances if KNOWN       
 
        function obj                 = VarianceModel(vr_type)
            obj.m_SetVarianceModel(vr_type);
        end

        % Set the error variance model (KNOWN/UNKNOWN)        
        function obj                 = m_SetVarianceModel(obj, vr_type)
            
            obj.d_VarianceModelType  = VarianceModelType.newType(vr_type);
            
            if (vr_type == "KNOWN")
               obj.d_VarianceComputationStatus = VarianceComputationStatusEnums.VARIANCE_KNOWN;
            else
                obj.d_VarianceComputationStatus = VarianceComputationStatusEnums.VARIANCE_TO_BE_COMPUTED;
            end
            
        end

        % Set the error variances (if KNOWN)
        function obj                    = m_setErrorVariances(obj, errorVariances)
            obj.d_errorVariances        = errorVariances;
            obj.d_defaultErrorVariances = errorVariances;
        end
        
        function obj                = m_resetVariancestoDefault(obj)
            obj.d_errorVariances    = obj.d_defaultErrorVariances;
        end

        % Get the error variances
        function errorVariances      = m_getErrorVariances(obj)
            errorVariances           = obj.d_errorVariances;
        end
        
%%  The functions in this block set the Variance - Covariance Mapper 
%   if the Error Variances are  UNKNOWN

        % Set the Variance Covariance Mapper if Error Variances are not known
        function obj                       = m_setMapper(obj, mapper_type)

%            if (obj.d_VarianceModelType == "VR_UnknownVariance")
                obj.d_varianceCovarianceMapper = VarianceCovarianceMapper(mapper_type);
%            else
%                disp("m_setMapper : Need not set Mapper if Variance Model is not UNKNOWN");
%            end
        end
        
        % Set Number of variables and lag only if its a DipcaMapper
        function obj = m_setNoOfVariablesAndLagInMapper(obj, no_of_variables, lag)
%            if (obj.d_varianceCovarianceMapper.m_getType() == "VCV_DipcaMapper")
                obj.d_varianceCovarianceMapper.m_setNoOfVariablesAndLag(no_of_variables, lag);
%            else
%                disp("m_setNoOfVariablesAndLagInMapper:The mapper type is not DIPCA");
%            end            
        end
        % set the no of constraints for dipca while estimating the order
        function obj = m_setNoOfConstraints(obj, no_of_constraints)
            obj.d_VarianceModelType.m_setNoOfConstraints(no_of_constraints);
        end
        
        % Set row indices and column indices only if its a GenericMapper
        function obj = m_setRowIndicesAndColumnIndices(obj, row_ind, col_ind)
            obj.d_varianceCovarianceMapper.m_setRowIndicesAndColumnIndices(row_ind, col_ind);
        end
        
        function covariance     = m_getCovarianceFromVariance(obj)
            covariance          = obj.d_varianceCovarianceMapper.m_mapVarianceToCovariance(obj.d_errorVariances);
        end
        
        
%%  Setting the Variance Solver when the Error Variances are UNKNOWN
%These functions to be swept into VR_Unknown??

        % Set a Variance Solver here
        function obj             = m_setVarianceSolver(obj, v_solver_type)
            disp("VarianceModel : m_setVarianceSolver ")
            obj.d_varianceSolver = VarianceSolver(v_solver_type);
        end
        
        % Compute the Variances here
        function obj             = m_solveForVariance(obj, data, no_of_constraints)
            disp("VarianceModel : m_solveForVariance ")
            obj.d_errorVariances
            [obj.d_errorVariances, obj.d_VarianceComputationStatus] = obj.d_varianceSolver.m_solveForVariance(data, obj.d_varianceCovarianceMapper, obj.d_errorVariances, no_of_constraints);
        end
        
%%  Computing the covariance (for a static case lag = 0)
        function [covariance, ret_code] = m_ComputeCovariance(obj, data, no_of_constraints )
           %Check if the d_VarianceModelType is UNKNOWN
           
            disp("VarianceModel :  m_ComputeCovariance")
            
            obj.m_resetVariancestoDefault();
            obj.m_solveForVariance(data, no_of_constraints);
            
            variance_vec        = obj.d_errorVariances;
            covariance          = obj.d_varianceCovarianceMapper.m_mapVarianceToCovariance(variance_vec);
            ret_code            = obj.d_VarianceComputationStatus;
        end
        
    end
    
end

