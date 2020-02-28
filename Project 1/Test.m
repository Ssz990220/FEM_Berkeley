k = 1;
L = 1;
E = 0.1;
syms x;
f = -14*x^2-k^2*cos(2*pi*k*x/L);
f = matlabFunction(f);
c_2 = 0.1+1/(4*pi^2*0.1);
c_1 = (1+10/(4*pi^2)-70/6-c_2);
u_true = -L^2/(4*pi^2*E)*cos(2*pi*k.*x/L)+7*x.^4/6/0.1+c_1.*x+c_2;
u_diff = diff(u_true);
u_diff = matlabFunction(u_diff);
u_true = matlabFunction(u_true);
figure;
x = 0:(32*100);
x = x*(1/(32*100));
curve_true = u_true(x);
plot(x, curve_true);
hold on;
n=0;
N_e = 4;
fem2 = FEM_1_D_Sizhe(L,N_e,E,2,'u-u',0.1,1,f,u_true,u_diff);
fem2.gen_K_R();
fem2.post_process();
fem2.solve();
fem2.gen_solution();
x = 0:(N_e*100);
x = x*(1/(N_e*100));
curve_true = u_true(x);
plot(x,fem2.solution);
hold on;
plot(fem2.Nodes,fem2.answer);
legend('u\_true','u\_post','u\_answer');