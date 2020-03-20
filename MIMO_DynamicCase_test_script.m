%% This is a test script for DIPCA for a MIMO system
% Input - data matrix (with variables along columns)
% Output - order, error variances, model, (confidence intervals calcaulation to be
% added)

clear
clc
%% Loading input and Initializing the XPCA object
% Input Data and Number of Inputs and Number of Outputs

%load 'test_files/test22.mat'

A = [0.67 0.67 0 0; 
    -0.67 0.67 0 0; 
    0 0 -0.67 -0.67; 
    0 0 0.67 -0.67];

B = [0.6598; 1.9698; 4.3171; -2.6436];
C = [-0.5749 1.0751 -0.5225 0.1830;   2.4027 0.7543 -0.2159 0.0982];
D = [-0.7139;  -0.5431];
E = [0; 0; 0; 0];
%E = [0.1762;0.5278;-0.5532;0.2983];

[b,a] = ss2tf(A,B,C,D);

%sys1=filt(b(1,:),a,1)
%sys2=filt(b(2,:),a,1)

%sys = [ sys1; sys2];
%sys = sys1;

N       = 1010; %Number of samples
n       = 4;    %State space order, used for data generation
m       = 2;   %output dimension
%m       = 1;
l       = 1;    %input dimension

% Input data
for k = 1:N
    f = 0;
    for j = 1:10
        f = f + sin(0.3898*pi*j*k);
    end
    u(:,k) = f; 
end

%Getting u and y matrices, x is state variable
x = zeros(n,N);
u = zeros(l,N);
y = zeros(m,N);

x(:,1) = 0;
for k = 1:N;
    f = 0;
    for j = 1:10;
        f = f + sin(0.3898*pi*j*k);
    end
    u(:,k) = f; 
    %+ 0.1*randn(1,1);
    if k==1
        x(:,k+1) = B*u(k);
        y(:,k) = D*u(k);
    else
        x(:,k+1) = A*x(:,k) + B*u(k) + E*randn(1);
        %x(:,k+1) = A*x(:,k) + B*u(k);
        y(:,k) = C*x(:,k)+D*u(k);
    end
end

% Output data
%y   = lsim( sys, u);
%f   = [ y u' ];

eu  = 0.1*randn(1,N);
%ey  = 0.2*randn(1,N);
ey  = [0.9*randn(1,N); 0.2*rand(1,N)];  

unew = u' + eu';
ynew = y' +  ey';
z = [ynew unew];
data = z;

g = xpca(data)

lag = 10;

%% The part handles the configuration of Model, Variance Model and Variance Solver

% Set the Model (Static/Dynamic)
g.m_setModel("DYNAMIC_MODEL")

g.m_setNoOfInputsAndOutputs(l,m)

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
order = g.m_computeOrder()

g.m_restackByOrder();
    
% Compute the Variance
% But variances are already computed in the previous step

g.m_computeCovarianceStackedByOrder();  
g.m_computeCovarianceStackedByLag();  

g.m_computeModel();

g.m_computeStateSpaceParameters()



