Nn = 64;
K = rand(Nn,Nn);
u0 = 5;
R = rand(Nn,1);
K(1,:) = 0;
K(:,1) = 0;
K(1,1) = 1;
u5 = 20;
K(5,:) = 0;
K(5,5) = 1;
R(1) = u0;
R(5) = u5;
a = K\R;

%%%2-D
K = zeros(Nn,Nn);
Ke = rand(4,4);
Conn = [1,12,24,48];
K(Conn,Conn) = Ke;


%%%vector_multiplication
a = rand(4,1);
b = rand(4,1);
c = a'*b;   %this method is faster than dot(a,b)


