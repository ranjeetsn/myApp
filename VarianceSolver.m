classdef VarianceSolver < handle
    % This class Solves for the variance given data and the solver type
    
    properties (Access = private)
        d_VarianceSolverType
    end
    
    methods
        
        function obj                = VarianceSolver(v_solver_type)
            obj.m_setVarianceSolver(v_solver_type);
        end
        
        function obj                =  m_setVarianceSolver(obj, v_solver_type)
           obj.d_VarianceSolverType =  VarianceSolverType.newType(v_solver_type);
        end
        
        function [variance, ret_code] = m_solveForVariance(obj, data, variance_covariance_mapper, initial_variance_vec, no_of_constraints)
            
            disp("In VarianceSolver : m_solveForVariance");
            
            [variance, ret_code]      = obj.d_VarianceSolverType.m_solveForVariance(data, variance_covariance_mapper, initial_variance_vec, no_of_constraints);
        end
        
    end
end

