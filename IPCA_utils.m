    classdef IPCA_utils
        % Class contains all the functions from IPCA 
        % Author : Shankar Narasimhan

        properties
        end

        methods (Static)

             function [x] = vectorize(M,rowind,colind)
                % converts the m x n matrix M into a column vector of length mn x 1 if
                % rowind/colind are not specified by arranging each column of M one below
                % the other starting with the first column.  Otherwise only the elements of M
                % corresponding to indices in rowind/colind are put in vector x.

                x = [];
                if ( nargin < 2 )
                    [m, n] = size(M);
                    for i = 1:n
                        x = [x; M(:,i)];
                    end
                else
                    if ( length(rowind) ~= length(colind) )
                        disp('Error: Length of rowind and colind are not equal');
                    end
                    for i = 1:length(rowind)
                        x = [x; M(rowind(i),colind(i))];
                    end
                end


             end

             function [M] = matricise(x, m, rowind, colind, flag)
                % Convert a vector x into a matrix.  If rowind and colind are not provided
                % then x is assumed to be of size mN elements which is converted to a form
                % of size m x N.  If rowind and colind are provided, then x is assumed to
                % same size as rowind/colind and these are assumed to be the nonzero
                % elements of a matrix of size m x m.  The location of the non-zeros are
                % assumed to be given by row and column indices in rowind/colind.  If flag
                % is one, then the matrix M is made symmetric by storing the same elements
                % in colind/rowind locations also.
                num_elements = length(x);
                if ( nargin < 5 )
                    flag = 0;   % default is M is not a symmetric matrix
                end
                if ( nargin < 3 )
                    %  Check if x contains mN elements
                    if ( mod(num_elements, m) )
                        disp('Error: Number of elements in vector X is not a multiple of m');
                    else
                        M = [];
                        N = floor(num_elements/m);
                        count = 1;
                        for i = 1:N
                            M = [M, x(count:count+m-1)];
                            count = count + m;
                        end
                    end
                else
                    if ( (num_elements ~= length(rowind)) | (num_elements ~= length(colind)) )
                        disp('Error: Size of x does not match with rowin/colind');
                    else
                        M = zeros(m,m);
                        for i = 1:num_elements
                            M(rowind(i),colind(i)) = x(i);
                        end
                        if ( flag )  % Symmetric matrix
                            for i = 1:num_elements
                                M(colind(i),rowind(i)) = x(i);
                            end
                        end

                    end
                end
            end


            function [Qe] = initialcove(A,Y,rowind,colind)
                %  function for generating initial estimate of error covariance matrix
                %  based on Keller's method (least squares solution)
                % INPUTS:
                % A : m x n constraint matrix
                % Y : data matrix n x N, n rows are variables and N columns are samples
                % rowind, colind : column vectors which contain the row and column indices corresponding to
                %                  non-zero elements of Qe.  The number of non-zero elements of Qe cannot
                %                  exceed m*(m+1)/2 for an identifiable problem. If rowind and colind
                %                  are not input Qe is assumed to be a diagonal matrix
                %
                % OUTPUTS
                %  Qe : Estimated measurement error covariance matrix
                %
                %  Check for inputs
                if nargin < 3
                    diagflag = 1;
                else
                    if ( size(rowind,1) )
                        diagflag = 0;
                    else
                        diagflag = 1;
                    end
                end
                [m n] = size(A);
                if ( ~diagflag )
                    nz = size(rowind,1);
                    maxnz = m*(m+1)/2;
                    if ( nz > maxnz )
                        disp ('The maximum number of nonzero elements of Qe that can be estimated exceeds limit');
                        return
                    end
                end
                % Construct D matrix
                G = [];
                for j = 1:n
                    C = [];
                    for i = 1:m
                        C = [C; A(i,j)*A(:,j)];
                    end
                    G = [G, C];
                end
                if ( ~diagflag )
                    for k = 1:nz
                        if ( rowind(k) ~= colind(k) )
                            C1 = [];
                            C2 = [];
                            for i = 1:m
                                C1 = [C1;A(i,rowind(k))*A(:,colind(k))];
                                C2 = [C2;A(i,colind(k))*A(:,rowind(k))];
                            end
                            G = [G, C1+C2];
                        end
                    end
                end
                % Construct RHS of equation
                nsamples = size(Y,2);
                R = A*Y;
                V = R*R'/nsamples;
                vecV = IPCA_utils.vectorize(V);
                % Least squares estimate
                estvar = abs((G'*G)\(G'*vecV));
                if ( diagflag )
                    Qe = diag(estvar);
                else
                    Qe = IPCA_utils.matricise(estvar,n,rowind,colind,1);
                end


            end

            %  function residuals for nonlinear equation solver

            function fres = obj_val(evar,A,Sr,rowind,colind)

                [ncon, nvar] = size(A);

                %  Calculate objective function value for every guess of the diagonal elements of
                %  square root of covariance matrix contained in vector evar

                fres = 0;
                if ( isempty(rowind) )
                    L = diag(evar);
                else
                    L = IPCA_utils.matricise(evar,nvar,rowind,colind,1);
                end

                % Compute objective function value (based on joint pdf of constraint
                % residuals
                Vmat = A*L*L'*A';
                Vinv = inv(Vmat);
                %     for j = 1:nsamples
                %         conres = A*Y(:,j);
                %         fres = fres + conres'*Vinv*conres;
                %     end
                fres = trace(Vinv*Sr) + log(det(Vmat));
            end

            
            function fres = obj_val1(evar,A,Sr, mapper)

                [ncon, nvar] = size(A);

                %  Calculate objective function value for every guess of the diagonal elements of
                %  square root of covariance matrix contained in vector evar

                fres = 0;
                
                L = mapper.m_mapVarianceToCovariance(evar);

                % Compute objective function value (based on joint pdf of constraint
                % residuals
                Vmat = A*L*L'*A';
                Vinv = inv(Vmat);
                fres = trace(Vinv*Sr) + log(det(Vmat));
            end
            
            
            function [u, sval, v, A, Qe, errflag] = ipca(Y, nfact, Agiven, rowind, colind, Qe)

                %                 IPCA.M   v. 2.0
                %
                % This function identifies a basis for the subspace orthogonal to Vx (space
                % in which true data vectors lie) given a sample of noisy measurements and
                % dimension of the subspace Vx.  The rows of matrix A form a basis for the
                % subspace orthogonal to Vx. It should be noted that the rows of A are not
                % orthonormal.  However the rows of A*L is orthonormal where L is the
                % cholesky factor of the error covarance matrix Qe.  Qe is also estimated
                % depending on the inputs specified by user.  The SVD of scaled data matrix
                % inv(L)*Y is also returned. Iterative PCA technique is used for this purpose.
                %
                % The measurement errors in each sample are assumed to have zero mean and
                % covariance matrix Qe. However, there could be correlation between
                % variables and across samples.  Furthermore, the true data vectors are
                % assumed to be deterministic (corresponds to functional model in
                % statistics). Unlike PCA which is optimal when the errors in all variables
                % are identical and are uncorrelated, the IPCA technique produces optimal
                % estimates for general covariance structures. 
                %
                % Function calls
                %
                % [u, sval, v, A, Qe, errflag] = ipca(Y, nfact) - diagonal Qe is estimated
                %
                % [u, sval, v, A, Qe, errflag] = ipca(Y, nfact, Agiven) - diagonal Qe is estimated,
                %                                First part of A is Agiven, remaining is estimated
                %
                % [u, sval, v, A, Qe, errflag] = ipca(Y, nfact, Agiven, rowind, colind) -
                % Qe is estimated.  Nonzero elements of Qe are given in rowind and colind
                %
                % [u, sval, v, A, Qe, errflag] = ipca(Y, nfact, Agiven, rowind, colind, Qe)
                % Qe is used for scaling data and is not estimated

                % INPUTS:
                % Y : data matrix n x N, n rows are variables and N columns are samples
                % nfact : dimension of subspace in which true data lies (number of
                %         independent variables).  A is of dimension (n - nfact) x N.  Theoretically
                %         if nfact is correct then there should be (n-nfact) singular values whose
                %         values are all close to unity.
                % Agiven : Rows of A which are already known (partial knowledge of
                %          constraint matrix A).  In this case the remaining rows of A are
                %          estimated. Specify empty matrix [] if no consttraints are known
                % rowind, colind : column vectors which contain the row and column indices corresponding to
                %                  non-zero elements of Qe.  The number of non-zero elements of Qe cannot
                %                  exceed (n-nfact)*(n-nfact+1)/2 for an identifiable problem.
                % Qe : if specified by user, then the IPCA technique uses it to identify A
                %      without iteration.  Otherwise Qe is estimated iteratively and returned
                % 
                % Notes:  rowind and colind can be set to empty sets [], if Qe is
                % specified.  If Qe is not specified, then Qe is assumed to be diagonal and
                % is estimated
                %
                % OUTPUTS
                % A : constraint matrix, rows of A form a basis for subspace orthogonal to
                %     true data space Vx.  First few rows of A will be Agiven and remaining are
                %     estimated
                % u : left singular vectors (orthonormal) of scaled data matrix 
                % sval : contains converged singular values of scaled data matrix
                % v : right singular vectors  (orthonormal) of scaled data matrix
                % Qe : estimated value of Qe if not specified by user.  Contains non-zero
                %      elements only in locations specified by user in rowind and colind.
                %      Otherwise only the diagonal elements are nozero. Otherwise same as Qe
                %      specified by user
                % errFlag indicates the termination conditions of the function;
                %             0 = normal termination (convergence)
                %             1 = maximum number of iterations exceeded
                %
                % The function can also produce a file - ipca.mat - if appropriate
                % lines are activated as indicated in the code.  This can be used to
                % follow convergence if desired.
                %
                % REFERENCES
                %
                %  1. Narasimhan, S. and S. L. Shah, "Model Identification and Error
                %  Covariance Matrix Estmation from Noisy Data using PCA," Control Engineering Practice. 16, 146-155 (2008).
                %
                %  2. Narasimhan, S. and N. Bhatt, "Deconstructing Prinicpal Component
                %  Analysis using a Data Recocniliation Perspective," Computers and
                %  Chemical Engineering, 77, 74-84 (2015).
                %
                %  3. P.D. Wentzell, D.T. Andrews, D.C. Hamilton, K. Faber, and
                %  B.R. Kowalski, "Maximum likelihood principal component analysis",
                %  J. Chemometrics, V 11, 339-366 (1997).
                %
                %  4. Fuller, W.A., Measurement Error Models, John Wiley & Sons, New York, 1987.
                %
                %-------------------------------------------------------------------------

                % Check for inputs and set defaults
                estflag = 0;
                diagflag = 0;
                conflag = 0;
                mgiven = 0;
                [nvar, nsamples] = size(Y);  % Number of measured variables & number of samples
                nsing = min(nvar,nsamples);
                if nargin > 2,
                    conflag = 1;
                    mgiven = size(Agiven,1);
                    if ( ~mgiven ) 
                        conflag = 0;
                    end
                end
                if nargin < 6,
                    estflag = 1;   % Qe needs to be estimated
                end
                if nargin < 4,
                    diagflag = 1;    % Qe structure is diagonal
                    rowind = [];
                    colind = [];
                end
                if ( estflag )
                    % Check if the number of non-zeros in Qe to be estimated satisfies
                    % limits
                    if ( diagflag )
                        nz = nvar;
                    else
                        nz = size(rowind,1);
                    end
                    maxnz = (nsing-nfact)*(nsing-nfact+1)/2;
                    if ( nz > maxnz )
                        disp ('The maximum number of nonzero elements of Qe that can be estimated exceeds limit');
                        errflag = 1;
                        return
                    end
                    if ( size(colind,1) )
                        if ( nz ~= size(colind,1) )
                            disp('The number of elements in rowind and colind do not agree');
                            errflag = 1;
                            return
                        end
                    end
                end

                sumsing = 0;  % Sum of singular values used to test for convergence
                sumold = 1;
                maxiter = 100;  % Maximum number of iterations
                vsmall = 1.0e-04; % Lower bound on elements of sqrt of Qe
                tol = 1.0e-04; % Relative tolerance for convergence (on sum of singular values)

                prop_factor = 0.01;  % Factor for initial estimate of error covariance matrix

                % Get an initial estimate of error covariance matrix from covariance matrix
                % of measurements
                if ( estflag )
                    SY = zeros(nvar);
                    for j = 1:nsamples
                        SY = SY + Y(:,j)*Y(:,j)';
                    end
                    meanY = mean(Y,2);
                    covY = (SY - nsamples*meanY*meanY')/(nsamples-1);
                % Set the initial estimates and bounds on error covariance matrix elements which
                % are estimated
                    if ( diagflag )
                        Qe = prop_factor*diag(diag(covY));
                        vlb = vsmall*ones(nvar,1);
                        vub = diag(SY);
                    else
                        for i = 1:nz
                            estvar(i)   = prop_factor*covY(rowind(i),colind(i));
                            vlb(i)      = vsmall;
                            vub(i)      = SY(rowind(i),colind(i));
                        end
                        Qe = IPCA_utils.matricise(estvar,nvar,rowind,colind,1);
                        %  Check if covE is singular
                        if ( rank(Qe) < nvar )
                            disp('Specified structure of Qe gives a singular matrix')
                            return
                        end
                    end
                end

                Nsqrt = sqrt(nsamples);
                if ( estflag )
                    L = diag(sqrt(diag(covY)));  % Autoscaling for PCA
                else
                    L = chol(Qe);
                    L = L';
                end
                Linv = inv(L);  % Store inverse of scaling matrix;
                            
                % Determine the residuals if part of the constraint matrix is given
                if ( conflag )
                    Yr = Y - Qe*Agiven'*inv(Agiven*Qe*Agiven')*Agiven*Y;
                else
                    Yr = Y;
                end
                                                    
                %  Determine initial estimate of A by scaling data using L
                Ys = Linv*Yr/Nsqrt;
                [u s v] = svd(Ys,0);
                % Get the matrix in terms of original variables
                Apart = u(:,nfact+1:nsing-mgiven)'*Linv;
                % Determine complete constraint matrix if part if specified
                if ( conflag )
                    A = [Agiven;Apart];
                else
                    A = Apart;
                end

                % Return if Qe is specified
                if ( ~estflag )
                    sval = diag(s);
                    errflag = 0;
                    return
                end

                % Iteratively estimate Qe and constraint matrix A
                iter = 0;
                % Use as initial estimate the cholesky factor of least squares estimate of error covariance
                % matrix by using Keller's solution

                if ( diagflag )

                    intial_x0 = IPCA_utils.initialcove(A,Y,rowind,colind);
                    x0 = sqrt(diag(intial_x0));

                else
                    x0 = IPCA_utils.vectorize(chol(IPCA_utils.initialcove(A,Y,rowind,colind))',rowind,colind);
                end

                while ( abs(sumsing-sumold) > sumold*tol )
                    iter = iter+1
                    %
                    % Compute residual covariance matrix and store for repeated use
                    %
                    res = A*Y;
                    Sr = cov(res',1);  % Normalize by N - number of samples    % Estimate maximum likelihood estimates of variables and error
                    % covariance matrix.
                    obj_fun = @IPCA_utils.obj_val;
                    [kest, fval] = fmincon(obj_fun,x0,[],[],[],[],vlb,vub,[],optimset('Display','iter','MaxFunEvals',50000),A,Sr,rowind,colind);
                    %     kest = [kest(1:40);ones(100,1)*kest(41);kest(42:length(kest))];
                             
                    % New estimates of the covariance matrix of errors
                    if ( diagflag )
                        L = diag(kest);
                    else
                        L = IPCA_utils.matricise(kest,nvar,rowind,colind,1);
                    end
                    Qe = L*L';
                    eig(A*Qe*A')
                    disp('Press any key to continue');
                    pause

                    % Determine the residuals if part of the constraint matrix is given
                    if ( conflag )
                        Yr = Y - Qe*Agiven'*inv(Agiven*Qe*Agiven')*Agiven*Y;
                    else
                        Yr = Y;
                    end
                    
                    %  Scale measurements using cholesky factor of new Qe
                    Linv = inv(L);
                    Ys = Linv*Yr/Nsqrt;

                    % Perform SVD on estimated variance matrix of X
                    [u s v] = svd(Ys,0);
                    % Get the matrix in terms of original variables
                    Apart = u(:,nfact+1:nsing-mgiven)'*Linv;
                    % Determine complete constraint matrix if part if specified
                    if ( conflag )
                        A = [Agiven;Apart];
                    else
                        A = Apart;
                    end

                    %  Determine the sum of singular values
                    sumold = sumsing
                    sval = diag(s)
                    sumsing = sum(sval(nfact+1:nsing))
                    if ( iter >= maxiter )
                        disp('Warning: Maximum iterations exceeded in IPCA');
                        sval = diag(s);
                        errflag = 1;
                        return
                    end
                             
                    % Set initial estimates of error covariances for next iteration to be
                    % the estimates from previous iteration
                    x0 = kest;
                end

                % Clean up and return
                sval = diag(s);
                errflag = 0;
            end


        end
    end


