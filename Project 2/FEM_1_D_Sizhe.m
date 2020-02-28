classdef FEM_1_D_Sizhe < handle

    properties
        L;      %total length
        N_n;    %number of nodes
        N_e;    %number of elements
        E;      %Young's Module
        Shape_function_type;
        Shape_func;
        B_C_type;
        B_C_1;
        B_C_2;
        f;      %
        K;      %Stiffness matrix
        R;      %Forcing Vector
        Nodes;  %Nodes' coordinate
        Conn;   %Connection of nodes
        answer; %weight of uN
        solution%fit curve
        u_true; %the true displacement
        u_diff;
    end
    
    methods
        function obj = FEM_1_D_Sizhe(L, N_e, E, Shape_function_type,B_C_type,B_C_1,B_C_2,f,u_true,u_diff)
            %FEM_1_D Construct an instance of this class
            %   Initialize elements, Stiff matrix and forcing vector
            obj.L = L;
            obj.N_e = N_e;
            obj.Shape_function_type = Shape_function_type;
            obj.N_n = Shape_function_type*N_e +1;
            obj.E = E;
            obj.Shape_func = Shape_function(Shape_function_type);
            obj.B_C_type = B_C_type;
            obj.B_C_1 = B_C_1;
            obj.B_C_2 = B_C_2;
            obj.f = f;
            obj.K = zeros(obj.N_n);
            obj.R = zeros(obj.N_n,1);
            obj.Nodes = zeros(obj.N_n,1);   %1-D
            obj.Conn = zeros(obj.N_e,2);    %1-D
            node_index = linspace(0,obj.N_e,obj.N_n);
            obj.Nodes = node_index'*obj.L/obj.N_e;
            obj.Conn(:,1)=node_index(1:obj.N_e)+1;
            obj.Conn(:,2)=node_index(2:obj.N_n)+1;
            obj.answer = [];
            obj.solution =[];
            obj.u_true = u_true;
            obj.u_diff = u_diff;
        end
        
        function obj = gen_K_R(obj)
            %Summary of this method goes here
            %   K is the stiffness matrix
            %   R is the forcing vector
            %   A is the solution
            for i = 1:obj.N_e
                id = obj.Conn(i,:);
                for l = 1:obj.Shape_func.num_of_Gaussian_Points
                    X = obj.Nodes(id)'*obj.Shape_func.shape_func(l,:)';
                    J = obj.Nodes(id)'*obj.Shape_func.shape_func_der(l,:)';
                    obj.K(id,id) = obj.K(id,id)+obj.Shape_func.wts(l)*obj.Shape_func.shape_func_der(l,:)'...
                        *1/J*obj.E*obj.Shape_func.shape_func_der(l,:)*1/J*J;
                    increase = obj.Shape_func.wts(l).*obj.f(X).*(obj.Shape_func.shape_func(l,:)'*J);
                    obj.R(id) = obj.R(id)+ increase;
                end
            end
        end
        
        function obj = post_process(obj)
            pass = true;
            while true
                switch obj.B_C_type
                case "u-t"
                    obj.K(1,:)=0;
                    obj.K(1,1)=1;
                    obj.R(1)=obj.B_C_1;
                    obj.R(end)=obj.R(end)+obj.B_C_2;
                case "u-u"
                    obj.K(1,:)=0;
                    obj.K(1,1)=1;
                    obj.K(end,:)=0;
                    obj.K(end,end)=1;
                    obj.R(1)=obj.B_C_1;
                    obj.R(end)=obj.B_C_2;
                case "t_u"
                    obj.K(end,:)=0;
                    obj.K(end,end)=1;
                    obj.R(end)=obj.B_C_2;
                    obj.R(1)=obj.R(1)+obj.B_C_1;
                otherwise
                    disp("Invalid Boundary Condition Type")
                    obj.reassign_BC();
                    pass = false; 
                end
                if pass
                    break
                end
            end
        end
        
        function obj = reassign_BC(obj,new_BC)
            obj.B_C_type = new_BC;
        end
        
        function obj = solve(obj)
            obj.answer = obj.K\obj.R;
        end
        
        function obj = gen_solution(obj)
           for i = 1:obj.N_e
               for j = 0:99
                   obj.solution((i-1)*100+j+1)=j/100*obj.answer(i+1)+(100-j)/100*obj.answer(i);
               end
               if i == obj.N_e
                   obj.solution(i*100+1)=obj.answer(i+1);
               end
           end
        end
        
        function value = get_value(obj,coor)
            length_of_element = obj.L/obj.N_e;
            id = fix(coor/length_of_element);
            zeta = (coor - (id+1/2)*length_of_element)/length_of_element;
            value = (1-zeta)/2*obj.answer(id)+(1+zeta)/2*obj.answer(id+1);
        end
        
        function error = get_error(obj)
            error_der = 0;
            error_num = 0;
            for i = 1:obj.N_e
                id = obj.Conn(i,:);
                for l = 1:obj.Shape_func.num_of_Gaussian_Points
                    X = obj.Nodes(id)'*obj.Shape_func.shape_func(l,:)';
                    J = obj.Nodes(id)'*obj.Shape_func.shape_func_der(l,:)';
                    dU_n = obj.answer(id)'*(obj.Shape_func.shape_func_der(l,:))'*(1/J);
                    error_der = error_der + obj.Shape_func.wts(l)*(obj.u_diff(X)-dU_n)*obj.E*...
                        (obj.u_diff(X)-dU_n)*J;
                    error_num = error_num + obj.Shape_func.wts(l)*(obj.u_diff(X))^2*obj.E*J;
                end
            end
            error = error_der/error_num;
        end
    end
end

