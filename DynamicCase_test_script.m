%% This is a test script for DIPCA for a SISO system
% Input - data matrix (with variables along columns)
% Output - order, error variances, model, (confidence intervals calcaulation to be
% added)

clear
clc
%% Loading input and Initializing the XPCA object
% Input Data and Number of Inputs and Number of Outputs

load 'dynamic_test_case.mat'

g = xpca(data);

lag = 15;

%% The part handles the configuration of Model, Variance Model and Variance Solver

% Set the Model (Static/Dynamic)
g.m_setModel("DYNAMIC_MODEL")

g.m_setNoOfInputsAndOutputs(1,1);

% Set the preprocessor
g.m_setPreprocessor("NULL")

% Preprocess the data
g.m_preprocessData()

    % Configure the Dynamic Model
    g.m_configureDynamicModel(lag)
   
% Stack the Data
g.m_stackDataByLag()

% Set the Variance Model (Known/ Unknown)
g.m_setVarianceModel("UNKNOWN")

    % If unknown, configure the settings for DIPCA Model 
    g.m_configureVarianceModelForDIPCA()

    % Or configure some genericmodel
    % g.m_configureVarianceModelForGenericMapper(row_index, column_index)

    g.m_setVarianceSolver("MLE_SOLVER")
    

% Set Order Calculator
g.m_initializeErrorVariances();
 
g.m_configureOrderCalculatorForHypothesisTesting()
%g.m_configureOrderCalculatorUserDefinedOrder(5)
  
%% This part handles the computation of order, error variances and the model

% Check for identifiability here??
order = g.m_computeOrder();

g.m_restackByOrder();
    
% Compute the Variance
% But variances are already computed in the previous step

g.m_computeCovarianceStackedByOrder();  
g.m_computeCovarianceStackedByLag();  
