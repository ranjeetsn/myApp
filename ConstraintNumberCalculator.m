classdef ConstraintNumberCalculator < handle
    %----------------------------------------------------------------------
    %%This class calculates the number of constraints in the given data
    % INPUT - data_obj (variables along the rows)
    %         variance_model
    %         d_min
    %         d_max
    % OUTPUT - Number of Constraints among the variables in the data matrix
    
    % ERROR - if returned d value is 9999999, then no "d" value in the range
    %         has been found
    %----------------------------------------------------------------------
    
    properties
        d_Data
        d_VarianceModel
        d_DMax
        d_DMin
        d_DValue
    end
    
    methods
%% These functions initialize the a ConstraintNumberCalculator object   

%         function obj        = ConstraintNumberCalculator(data, variance_model, d_min, d_max)
%             disp('In ConstraintNumberCalculator: Constructor Overloaded');
%             obj.m_setData(data);
%             obj.m_setVarianceModel(variance_model);
%             obj.m_setDMin(d_min);
%             obj.m_setDMax(d_max);
%         end
%         
        function obj        = m_setData(obj, data)
            obj.d_Data      = data;
        end
        
        function obj        = m_setVarianceModel(obj, variance_model)
            obj.d_VarianceModel = variance_model;
        end
        
        function obj        = m_setDMax(obj, d_max)
            obj.d_DMax      = d_max;
        end
        
        function obj        = m_setDMin(obj, d_min)
            obj.d_DMin      = d_min;
        end
%% These functions handle the computation of the number of constraints

        function [number_of_constraints, eig_values, ret_code] = m_ComputeNumberOfConstraints(obj)
            
            noise_model         = NoiseModel.getInstance();
            Z                   = obj.d_Data.m_getData();
            no_of_realizations  = size(Z, 2);
            
            if( noise_model.m_getNoiseModel() == NoiseModelEnums.NOISELESS)
                
                % Compute number of constraints for noiseless case
                eig_values = SysidUtils.m_computeEigenValues(Z);
                number_of_constraints     = SysidUtils.HypothesisTest1(eig_values, no_of_realizations, SysidUtils.ALPHA_HYPOTHESIS_TESTING());
                obj.d_DValue              = number_of_constraints;
                ret_code                  = ReturnCodeEnums.SUCCESS;
                return;
                
            else
                
                for d_val = obj.d_DMax : -1 : obj.d_DMin
                    
                    [eig_values, rcode] = SysidUtils.m_ComputeEigenvaluesForGivenD( obj.d_Data,  d_val, obj.d_VarianceModel);
                    
                    if (rcode ~= VarianceComputationStatusEnums.VARIANCE_COMPUTATION_FAILED_DUE_TO_UNKNOWN_REASON && ...
                            rcode ~= VarianceComputationStatusEnums.VARIANCE_COMPUTATION_FAILED_DUE_TO_NON_IDENTIFIABILITY)
                        
                        if ( SysidUtils.HypothesisTestForaGivenD(eig_values, no_of_realizations, SysidUtils.ALPHA_HYPOTHESIS_TESTING(), d_val ) )
                            
                            number_of_constraints   = d_val;
                            obj.d_DValue            = d_val;
                            ret_code                = ReturnCodeEnums.SUCCESS;
                            return;
                        end
                        
                    end
                    
                end
            end
            number_of_constraints   = obj.DEFAULT_NUMBER_OF_CONSTRAINTS(); 
            obj.d_DValue            = obj.DEFAULT_NUMBER_OF_CONSTRAINTS();
            eig_values              = [];
            ret_code                = ReturnCodeEnums.FAILED;
        end
        
    end
    methods (Static)
        %% Default value of number of constraints when the computation fails      

        function d          = DEFAULT_NUMBER_OF_CONSTRAINTS(obj)
            d               = 9999999;
        end
        
    end
    
 %   methods (Access = private)


%% Helper functions to compute variances and eigen values     
% 
%         function eig_vals   = m_ComputeEigenvaluesForGivenD(obj, d_val)
%             disp("OC_HypothesisTesting: m_ComputeEigenvaluesForGivenD")
%             co_variance     = obj.m_ComputeNoiseCovarianceForGivenD(d_val);
%             eig_vals        = obj.m_ComputeEigenValuesGivenDataAndCovariance(co_variance);
%             
%         end
%         
%         function covariance = m_ComputeNoiseCovarianceForGivenD(obj, d_val)
%             obj.d_VarianceModel.m_setNoOfConstraints(d_val);
%             covariance      = obj.d_VarianceModel.m_ComputeCovariance(obj.d_Data, d_val);
%         end
%         
%         function eig_vals   = m_ComputeEigenValuesGivenDataAndCovariance(obj,  co_variance)
%             eig_vals        = SysidUtils.m_ComputeEigenValuesforGivenDataAndCovariance(obj.d_Data.m_getData(), co_variance);
%         end
%     end
%     
end

