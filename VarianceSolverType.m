classdef VarianceSolverType < handle
    % This is an interface/ abstract class for variance solver classes
    % Any class that is implementing a variance solver needs to inherit
    % this class
    
    properties
    end
    
    
    methods (Abstract)
        [variance, ret_code] = m_solveForVariance(obj, data, variance_covariance_mapper);
    end
    
    
    methods (Static)
        
        function type = newType(v_solver_type)
            
            switch v_solver_type
                
                case 'TLS_SOLVER'
                    type = VS_TLSSolver;
                    
                case 'MLE_SOLVER'
                    type = VS_MLESolver;
                    
                otherwise
                    error('Not a valid Variance Solver Type');   
            end
            
        end
        
    end

end

