classdef VCV_GenericMapper < VarianceCovarianceMapperType
    % This is a concerte class implementing a Variance Covariance Mapper 

    properties (Access = private)
        d_rowInd = [];
        d_colInd = [];
    end
    
    methods
        
        function m_setRowIndexAndColumnIndex(obj, row_ind, col_ind)
            % Basic sanity check to see if the size of both 
            % row_ind and col_ind are the same
            disp(" In VCV_GenericMapper : m_setRowIndexAndColumnIndex");
            
            % Check row index is a column vector
            if (size(row_ind,2) ~= 1)
                disp("row_ind is not a column vector: check the dimension");
                return;
            end
            
            % Check if the sizes of row and column index are same before
            % assigning to the object properties
            if ( (size(row_ind,1) == size(col_ind, 1)) && ...
                 (size(row_ind,2) == size(col_ind, 2))  ) 
                  obj.d_rowInd = row_ind;
                  obj.d_colInd = col_ind;
            else
                  disp("row and column index dimensions are not equal");
            end
            
        end
        
        
        % Maps the variance vector to the diagonal matrix of covariance
        
        function covariance_matrix = Map(obj, variance_vector)
          % Check to see if the size of the variance_vector
          % is same as that of row_ind 
          
          disp(' In the  VCV_GenericMapper - Map routine')
          
          
          if (  (size(obj.d_rowInd, 1) == size(variance_vector, 1)) && ...
                (size(obj.d_rowInd, 2) == size(variance_vector, 2))  ) 
                % start building the covariance matrix
                
                No_of_Elements = size(obj.d_rowInd, 1);
                
                for i = 1:No_of_Elements
                    
                    row_index    = obj.d_rowInd(i, 1);
                    column_index = obj.d_colInd(i, 1);
                    
                   covariance_matrix(row_index, column_index) =  variance_vector(i, 1);
                end
                
           else
                disp("row index and variance_vector  dimensions are not equal");
           end
          
          

        end
        
    end
    
end

