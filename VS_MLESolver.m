classdef VS_MLESolver < VarianceSolverType
    % This is a concerte class implementing a Variance Solver Type using
    % MLE formulation
    properties
    end
    
    methods
        
        % Solves for variance given data matrix and the mapper
        
        function [variance, ret_code] = m_solveForVariance(obj, data, variance_covariance_mapper, initial_variance_vec, no_of_constraints)
            
          disp(' In VS_MLESolver : m_solveForVariance')
          variance      = [];
          
          nvar          = size(initial_variance_vec, 1);
          Y             = data.m_getData();
          no_of_samples = size(Y,2);
          Nsqrt         = sqrt(no_of_samples);
          sum_sing_new  = 0;
          sum_sing_old  = 1;
%           maxiter       = 100;
          vsmall        = 1.0e-04;
          tol           = 1.0e-04;
          iter          = 0;

          vlb           = vsmall * ones(nvar,1);
          No_of_Unknowns = size(initial_variance_vec, 1);
          
          % Call the identifiability check module
          if ( (no_of_constraints * (no_of_constraints + 1) / 2) < No_of_Unknowns)
              ret_code = VarianceComputationStatusEnums.VARIANCE_COMPUTATION_FAILED_DUE_TO_NON_IDENTIFIABILITY;
              return;
          end
          
          x0    = initial_variance_vec;
          L     = variance_covariance_mapper.m_mapVarianceToCovariance(x0); 
          Qe    = L * L';
          
          while ( abs(sum_sing_new - sum_sing_old) > tol)
          
            iter    = iter + 1;
            
            Ys      = SysidUtils.ConstructScaledZ1(Y, Qe);
            Ys      = Ys/Nsqrt;
            [u s v] = SysidUtils.GetSVD(Ys);
            Ak      = SysidUtils.GetConstraintMatrix(u, no_of_constraints);
            Linv    = inv(chol(Qe,'lower'));
            A       = Ak' * Linv ;
           
            res     =  A * Y ;
            Sr      = cov(res' , 1);
              
            obj_fun                             = @IPCA_utils.obj_val1;
            [kest, fval, exitflag, output]      = fmincon(obj_fun,x0,[],[],[],[],vlb,[],[],                 ...
                                                            optimset('Display','iter','MaxFunEvals',50000), ...
                                                            A,Sr,variance_covariance_mapper);
           
            L  = variance_covariance_mapper.m_mapVarianceToCovariance(kest); %Function name needs to change. Mapping to Cholesky.            
            
            Qe = L * L'  
            
            eig( A * Qe * A' )
            
            if ( iter > XpcaConstants.MAX_VARIANCE_ITERATIONS)
                disp (' maximum iterations hit: exiting variance computation')
                %ret_code =
                %VarianceComputationStatusEnums.VARIANCE_COMPUTATION_MAX_ITERATIONS_HIT;
                ret_code = VarianceComputationStatusEnums.VARIANCE_COMPUTATION_FAILED_DUE_TO_UNKNOWN_REASON;
                return;
            end

            sum_sing_old    = sum_sing_new;   
            no_of_variables = size(Qe, 1);
            sval            = diag(s);
            sum_sing_new    = sum(sval(no_of_variables - no_of_constraints + 1:no_of_variables));
            x0              = kest;
            
          end
%          A
          variance = x0.^2;  %This might need to be changed later
          ret_code = VarianceComputationStatusEnums.VARIANCE_COMPUTED_SUCCESSFUL;
        end
        
    end
    
end

