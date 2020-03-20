%% This script is used to test the class to calculate the no of constraints
% INPUT     - data (variables along columns), variance model
% OUTPUT    - number of constraints
%           - eigen values

clear
clc

% Setting up the Data object
load InputdataIPCA.mat
data        =   Z;
data_obj    =   Dataobj(data'); % The variables are rows and realizations are columns

% Setting up the Variance Model
v           =   VarianceModel("UNKNOWN");
v.m_setVarianceSolver("MLE_SOLVER")
v.m_setMapper("DIAGONAL_MAPPER")

nvar        =   size(data_obj.m_getData(), 1);
var         =   ones(nvar, 1);
v.m_setErrorVariances(var)

% Set limits for d_value
d_max       =   5;
d_min       =   2;

%struct(v)
%struct(d)

% Initialize the ConstraintNumberCalculator and compute 
%cnc                                 = ConstraintNumberCalculator(data_obj, v, d_min, d_max);
cnc                                 = ConstraintNumberCalculator();
cnc.m_setData(data_obj);
cnc.m_setVarianceModel(v);
cnc.m_setDMin(d_min);
cnc.m_setDMax(d_max);

[no_of_constraints, eigen_values]   = cnc.m_ComputeNumberOfConstraints()
