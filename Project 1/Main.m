k=1;
syms x;
L = 1;                %Total length
N_e = 256;              %number of elements
E = 0.1;
f = -14*x^2-k^2*cos(2*pi*k*x/L);
f = matlabFunction(f);
c_2 = 0.1+1/(4*pi^2*0.1);
c_1 = (1+10/(4*pi^2)-70/6-c_2);
u_true = -L^2/(4*pi^2*E)*cos(2*pi*k.*x/L)+7*x.^4/6/0.1+c_1.*x+c_2;
u_true = matlabFunction(u_true);

fem = FEM_1_D(1,N_e,0.1,1,'u-u',0.1,1,f);
fem.gen_K_R();
fem.post_process();
fem.solve();
fem.gen_solution();

x = 0:(N_e*100);
x = x*(1/(N_e*100));
curve_true = u_true(x);
figure;
plot(x,curve_true,x,fem.solution);
legend("U\_True","FEM");

