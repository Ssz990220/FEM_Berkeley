classdef FEM_1_D < handle
    %FEM_1_D is a class which solves 1-D FEM problem
    
    properties
        L;      %total length
        N_n;    %number of nodes
        N_e;    %number of elements
        E;      %Young's Module
        Shape_function_type;
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
    end
    
    methods
        function obj = FEM_1_D(L, N_e, E, Shape_function_type,B_C_type,B_C_1,B_C_2,f)
            %FEM_1_D Construct an instance of this class
            %   Initialize elements, Stiff matrix and forcing vector
            obj.L = L;
            obj.N_n = N_e + 1;
            obj.E = E;
            obj.N_e = N_e;
            obj.Shape_function_type = Shape_function_type;
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
            obj.Conn(:,1)=node_index(1:obj.N_e);
            obj.Conn(:,2)=node_index(2:obj.N_n);
            obj.answer = [];
            obj.solution =[];
        end
            
            
        function obj = gen_K_R(obj)
            %Summary of this method goes here
            %   K is the stiffness matrix
            %   R is the forcing vector
            %   A is the solution
            for i = 1:obj.N_e
                %Local Stiffness Metrix
                length_of_element = obj.Nodes(i+1)-obj.Nodes(i);
                Jacob = length_of_element / 2;
                k_local_11 = (-1/2)*obj.E*(-1/2)*1/Jacob*2;
                k_local_12 = (-1/2)*obj.E*(1/2)*1/Jacob*2;
                k_local_21 = (1/2)*obj.E*(-1/2)*1/Jacob*2;
                k_local_22 = (1/2)*obj.E*(1/2)*1/Jacob*2;
                obj.K(i,i) = obj.K(i,i)+k_local_11;
                obj.K(i,i+1) = obj.K(i,i+1)+k_local_12;
                obj.K(i+1,i) = obj.K(i+1,i)+k_local_21;
                obj.K(i+1,i+1) = obj.K(i+1,i+1)+k_local_22;
                zeta = 0.577350269189626;
                x_local_1 = obj.Nodes(i)*(1-zeta)/2+obj.Nodes(i+1)*(1+zeta)/2;
                x_local_2 = obj.Nodes(i)*(1+zeta)/2+obj.Nodes(i+1)*(1-zeta)/2;
                r_local_1 = obj.f(x_local_1)*(1-zeta)/2*Jacob+obj.f(x_local_2)*(1-zeta)/2*Jacob;
                r_local_2 = obj.f(x_local_1)*(1+zeta)/2*Jacob+obj.f(x_local_2)*(1+zeta)/2*Jacob;
                obj.R(i) = obj.R(i)+r_local_1;
                obj.R(i+1) = obj.R(i+1) + r_local_2;  
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
    end
end

