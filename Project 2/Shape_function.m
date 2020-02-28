classdef Shape_function < handle
    %SHAPE_FUNCTION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        shape_func;
        shape_func_der;
        wts;
        num_of_points;
        Gaussian_points;
        num_of_Gaussian_Points;
    end
    
    methods
        function obj = Shape_function(order_of_polynomial)
            obj.num_of_points = order_of_polynomial;
            obj.num_of_Gaussian_Points = ceil((order_of_polynomial+1)/2)+1;
            switch obj.num_of_Gaussian_Points
                case 2
                    obj.wts = [1,1];
                    obj.Gaussian_points = [0.577350269189626;-0.577350269189626];
                case 3
                    obj.wts = [0.888888888888889,0.555555555555556,0.555555555555556];
                    obj.Gaussian_points = [0;0.774596669224148;-0.774596669224148];
                case 4
                    obj.wts = [0.652145154862546,0.347854845137454,0.652145154862546,0.347854845137454];
                    obj.Gaussian_points = [0.339981043584856;0.861136311594053;...
                        -0.339981043584856;-0.861136311594053];
                case 5
                    obj.wts = [0.568888888888889;0.478628670499366;0.236926885056189;...
                        0.478628670499366;0.236926885056189];
                    obj.Gaussian_points = [0;0.538469310105683;0.906179845938664;...
                        -0.538469310105683;-0.906179845938664];
            end
            switch order_of_polynomial
                case 1
                    obj.shape_func = [(1-obj.Gaussian_points)./2,(1+obj.Gaussian_points)./2];
                    obj.shape_func_der = [-1/2,1/2;-1/2,1/2];
%                 case 2
%                     pass
            end
        end
    end
end

