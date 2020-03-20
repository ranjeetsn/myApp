classdef Dataobj
    % Input data matrix is stored in this object
    
    properties (Access = private)
        d_data
    end
    
    methods

        % Initialize the Data object
        function obj = Dataobj(input_data)
            obj.d_data          = input_data;
        end

        % Get the data
        function data = m_getData(obj)
                data = obj.d_data;
        end

    end
end

