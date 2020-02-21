L = 1;                %Total length
E = 0.1;
shape_func = Shape_function(1);
error =zeros(2,5,6);

for N = 0:5
    k = 2^N;
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
    for n = 0:4
        N_e = k*4*2^n;
        fem2 = FEM_1_D_Roger(L,N_e,E,1,'u-u',0.1,1,f,u_true,u_diff);
        fem2.gen_K_R();
        fem2.post_process();
        fem2.solve();
        fem2.gen_solution();
        error(1,n+1,N+1) = N_e;
        error(2,n+1,N+1) = fem2.get_error();
        x = 0:(N_e*100);
        x = x*(1/(N_e*100));
        curve_true = u_true(x);
        plot(x,fem2.solution);
        hold on;
    end
    title(append('k=',string(k)));
    legend("U\_True",...
    append("N_e = ",string(k*4*2^0)," Error = ",string(error(2,1,N+1))),...
    append("N_e = ",string(k*4*2^1)," Error = ",string(error(2,2,N+1))),...
    append("N_e = ",string(k*4*2^2)," Error = ",string(error(2,3,N+1))),...
    append("N_e = ",string(k*4*2^3)," Error = ",string(error(2,4,N+1))),...
    append("N_e = ",string(k*4*2^4)," Error = ",string(error(2,5,N+1))));
    hold off;
end

error_log = log(error);
error_log(1,:,:) = -error_log(1,:,:);
fit = zeros(6,2);
figure;

for i = 1:6
    f = polyfit(error_log(1,:,i),error_log(2,:,i),1);
    plot(error_log(1,:,i),error_log(2,:,i));
    fit(i,:) = f;
    hold on;
end

title('error to 1/N_e');
legend(append('k=1, y=',string(fit(1,1)),'+',string(fit(1,2))),...
    append('k=2, y=',string(fit(2,1)),'+',string(fit(2,2))),...
    append('k=4, y=',string(fit(3,1)),'+',string(fit(3,2))),...
    append('k=8, y=',string(fit(4,1)),'+',string(fit(4,2))),...
    append('k=16, y=',string(fit(5,1)),'+',string(fit(5,2))),...
    append('k=32, y=',string(fit(6,1)),'+',string(fit(6,2))));
xlabel('log(1/N_e)');
ylabel('log(error)');
