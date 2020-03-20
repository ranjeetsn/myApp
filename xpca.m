classdef xpca < handle
    
% This class handles various versions of PCA (PCA, PCA on scaled data, IPCA, DPCA, DIPCA)
% for determining the model, (order, error variances) 

% Multiple ways to preprocess the input data
% Handle static and dynamic models
% Error variances known/unknown or homoskedastic/heteroskedastic
% Calculate order of a dynamic system


    properties (Access = private)

        % Interfaces
        d_Preprocessor
%        d_NoiseModel  % Noise Model has been defined as Singleton for
                       % Universal Access
        d_Model 
        d_OrderCalculator
        d_VarianceModel
        d_ConfidenceIntervalCalculator
        
        % Data
        d_Dataobj  
        d_PreprocessedData
        d_StackedByLag
        d_StackedByOrder
        d_Covariance
               
    end
    
    
    methods
%% Creating the basic object by setting the data, number of inputs and outputs        
        % Initialize the xpca object with the input data matrix
        function obj                = xpca(data)
            obj.d_Dataobj           = Dataobj(data);
        end


%% Setting Preprocessor and Model        
        % Set the preprocessor using this method 
        function obj                = m_setPreprocessor(obj, pp_type)
            obj.d_Preprocessor      = Preprocessor(pp_type);
        end
        
        % Set the model(Static/Dynamic) using the method below
        function obj                = m_setModel(obj, model_type)
            obj.d_Model             = Model(model_type);
            obj.d_Model.m_setData(obj.d_Dataobj);
        end
%% Setting the Noise Modl
    function obj                    = m_setNoiseModel(obj, noise_model_type)
        noise_model = NoiseModel.getInstance();
        noise_model.m_setNoiseModel(noise_model_type);
    end
        
%% If the model is DYNAMIC configure some settings
        function obj                = m_configureDynamicModel(obj, lag)
            obj.d_Model.m_setLag(lag);
        end
        
        % Set the number of inputs and number of outputs in the Data object
        function obj                = m_setNoOfInputsAndOutputs(obj, noOfInputs, noOfOutputs)
            obj.d_Model.m_setNoOfInputsAndOutputs(noOfInputs, noOfOutputs);
        end
        
%% Setting the Variances if KNOWN or configure relevant settings in the next block if UNKNOWN        
        % Set the error variances here if known or unknown (compute using
        % iterative method)
        function obj                = m_setVarianceModel(obj, variance_model)
            obj.d_VarianceModel     = VarianceModel(variance_model);
        end
        
        function var_model          = m_getVarianceModel(obj)
            var_model               = obj.d_VarianceModel;
        end

        % Set the error variances if known
        function obj                = m_setErrorVariances(obj, var_array)
            obj.d_VarianceModel.m_setErrorVariances(var_array);
        end
        
%% The below block is used to configure the settings of the mapper and solver used to estimate the covariance matrix        

%% Configure the Mapper in this block

        % Configure the Diagonal Mapper for the Static case
        function obj                = m_configureVarianceModelForDiagonalMapper(obj)
            obj.d_VarianceModel.m_setMapper("DIAGONAL_MAPPER");
        end

        % Configure the Variance-Covariance Mapper for DIPCA 
        function obj                = m_configureVarianceModelForDIPCA(obj)
            
            obj.d_VarianceModel.m_setMapper("DIPCA_MAPPER");
            [n_inputs, n_outputs]   = obj.d_Model.m_getNoOfInputsAndOutputs();
            no_of_variables         = n_inputs + n_outputs;
            lag                     = obj.d_Model.m_getLag();
            obj.d_VarianceModel.m_setNoOfVariablesAndLagInMapper(no_of_variables, lag);
%            obj.d_Model.m_setMapper(obj.d_VarianceModel.d_varianceCovarianceMapper);
        end
        
        % Configure Generic Mapper by row ind and col ind
        function obj                = m_configureVarianceModelForGenericMapper(obj, row_ind, col_ind)
            obj.d_VarianceModel.m_setMapper("GENERIC_MAPPER");
            obj.d_VarianceModel.m_setRowIndicesAndColumnIndices(row_ind, col_ind);
        end
        
%% Configure the Variance solver in this block

        % Set the Variance Solver here
        function obj                = m_setVarianceSolver(obj, v_solver_type)
             obj.d_VarianceModel.m_setVarianceSolver(v_solver_type);
        end
        
%% Set the Order Calculator in this block
    
        function obj                = m_configureOrderCalculatorUserDefinedOrder(obj, order)
            obj.d_OrderCalculator   = OrderCalculator("USER_DEFINED");
            obj.d_OrderCalculator.m_configureOrderCalculator( obj.d_Dataobj, ...
                                                              obj.d_Model, ...
                                                              obj.d_VarianceModel);
            obj.d_OrderCalculator.m_setOrder(order);

        end

        function obj                = m_configureOrderCalculatorForHypothesisTesting(obj)
            obj.d_OrderCalculator   = OrderCalculator("HYPOTHESIS_TESTING");
            obj.d_OrderCalculator.m_configureOrderCalculator( obj.d_Dataobj, ... 
                                                              obj.d_Model, ...
                                                              obj.d_VarianceModel);

        end

%% With all the settings and configuration done in the blocks above, the below block will handle the computation

        % Preprocess data here  
        function obj                = m_preprocessData(obj)
             obj.d_PreprocessedData = Dataobj(obj.d_Preprocessor.PreprocessData(obj.d_Dataobj.m_getData()));
        end
         
        % Stack the data here 
        function obj                = m_stackDataByLag(obj)
             obj.d_StackedByLag     = Dataobj(obj.d_Model.m_StackDataByLag(obj.d_PreprocessedData));
        end
        
        function [d_order, eigen_values]    = m_computeOrder(obj)
            [d_order, eigen_values]         = obj.d_OrderCalculator.m_computeOrder();
            obj.d_Model.m_setOrder(d_order);
        end
        
        function obj                = m_restackByOrder(obj)
            order                   = obj.d_OrderCalculator.m_getOrder();
            obj.d_StackedByOrder    = Dataobj(obj.d_Model.m_StackDataByOrder(obj.d_PreprocessedData));            
        end
         
        % Compute Error Variances here if not known
        % Initialize error variances to zeros - can be changed to 
        % initialize to keller's solution
        
        function obj                = m_initializeErrorVariances(obj)
            
            no_of_vars              = obj.d_Model.m_getNoOfVariables();
            var_array               = ones(no_of_vars,1);           
            obj.m_setErrorVariances(var_array);
        end
        
        function obj                = m_computeCovarianceStackedByLag(obj)
            
            lag                     = obj.d_Model.m_getLag();
            order                   = obj.d_Model.m_getOrder();
            no_of_constraints       = lag - order + 1;
            obj.d_VarianceModel.d_varianceCovarianceMapper.m_setLag(lag);
            obj.d_Covariance        = obj.d_VarianceModel.m_ComputeCovariance(obj.d_StackedByLag, no_of_constraints);
        end
        
        function obj                = m_computeCovarianceStackedByOrder(obj)
            
            [no_of_inputs, no_of_outputs]   = obj.d_Model.m_getNoOfInputsAndOutputs();
            order                   = obj.d_Model.m_getOrder();
            no_of_constraints       = (no_of_outputs - 1) * order + no_of_outputs;
            obj.d_VarianceModel.d_varianceCovarianceMapper.m_setLag(order);
            obj.d_Covariance        = obj.d_VarianceModel.m_ComputeCovariance(obj.d_StackedByOrder, no_of_constraints);
        end
        
        function [Covariance, eig_vals] = m_computeCovariancesForStaticCase(obj)
            
            data                    = obj.d_Dataobj.m_getData();
            cnc                     = SysidUtils.m_ConfigureConstraintNumberCalculatorForStatic(Dataobj(data'), obj.d_VarianceModel);
            [num_equal_eigs, eig_vals] = cnc.m_ComputeNumberOfConstraints();

            no_of_constraints       = num_equal_eigs;
            obj.d_Covariance        = obj.d_VarianceModel.m_ComputeCovariance(Dataobj(data'), no_of_constraints)
            Covariance              = obj.d_Covariance
        end
        
        % Compute the Model using the data matrix and error variances
        function model              = m_computeModel(obj)
            
           noise_model              = NoiseModel.getInstance();
           
           if (noise_model.m_getNoiseModel() == NoiseModelEnums.NOISY)
               variance_vec         = obj.d_VarianceModel.m_getErrorVariances();
               variance_mapper      = obj.d_VarianceModel.d_varianceCovarianceMapper;
               % scale the data matrix with the inverse of cholesky factor and perform SVD
               model = obj.d_Model.m_computeModel(obj.d_Dataobj, variance_vec, variance_mapper);
           else
               model = obj.d_Model.m_computeModel(obj.d_Dataobj);
           end
           
        end
        
        function confidence_intervals = m_computeConfidenceIntervals(obj, number_of_bootstraps, size_of_bootstraps, alpha)
            
            cic                       = SysidUtils.m_ConfigureConfidenceIntervalCalculator(obj.d_StackedByOrder, number_of_bootstraps, ...
                                                                                            size_of_bootstraps, alpha);
            
            order                     = obj.d_Model.m_getOrder();
            [n_inputs, n_outputs]     = obj.d_Model.m_getNoOfInputsAndOutputs();
            
            noise_model = NoiseModel.getInstance();
            if (noise_model.m_getNoiseModel() == NoiseModelEnums.NOISY)
                variance_vec              = obj.d_VarianceModel.m_getErrorVariances();
                variance_mapper           = obj.d_VarianceModel.d_varianceCovarianceMapper;   
                confidence_intervals      = cic.m_computeConfidenceIntervals(order, n_outputs, order, variance_vec, variance_mapper);
            else
                 confidence_intervals      = cic.m_computeConfidenceIntervals(order, n_outputs, order);
            end
        end
        
        function m_computeStateSpaceParameters(obj)
            obj.d_Model.m_computeStateSpaceParameters();
        end
         
    end
    
end
