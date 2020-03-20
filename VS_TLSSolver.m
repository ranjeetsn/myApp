classdef VS_TLSSolver < VarianceSolverType
    % This is a concerte class implementing a Variance Solver Type using
    % TLS method
    properties
    end
    
    methods
        
        % Solves for variance given data matrix and the mapper
        
        function [variance, ret_code] = m_solveForVariance(obj, data, variance_covariance_mapper, initial_variance_vec, no_of_constraints)
            
          disp(' In VS_TLSSolver : m_solveForVariance')
          % Keller's solution here??
          variance = [];

        end
        
    end
    
end

