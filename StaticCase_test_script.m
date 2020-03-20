%% This file is the script to set up and test the Static case
%
% Input     - data file
% Output    - covariance & model
%
%% Loading the data and Initializing the XPCA object 

load InputdataIPCA.mat

data = Z;

g = xpca(data)

%% This part handles the configuration of Model, Variance Model, Variance Solver

g.m_setModel("STATIC_MODEL")

g.m_setPreprocessor("NULL")

g.m_preprocessData()

g.m_setVarianceModel("UNKNOWN")

g.m_configureVarianceModelForDiagonalMapper()

g.m_setVarianceSolver("MLE_SOLVER")

g.m_configureOrderCalculatorForHypothesisTesting(

% Initial guess of variance to start with
g.m_initializeErrorVariances();

%% The actual computation of the variances and model

covariance  = g.m_computeCovariancesForStaticCase()

model       = g.m_computeModel()