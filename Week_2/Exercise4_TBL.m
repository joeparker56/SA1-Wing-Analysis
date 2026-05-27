clear;
close all;


global re_l dudx u0

re_l = 1e7;
dudx =  -0;
u0 = 1;


x0 = 0.01;
thick0(1) = 0.023*x0*(re_l*x0)^(-1/6);
thick0(2) = 1.83*thick0(1);

[delx, thickhist] =ode45(@thickdash, [0, 0.99], thick0);

theta7 = 0.037*delx.*(re_l.*(delx+x0)).^(-1/5);
theta9 = 0.023*delx.*(re_l.*(delx+x0)).^(-1/6);

figure;
plot(delx+x0, thickhist(:,1))
hold on;
plot(delx+x0, theta7)
hold off;
hold on;
plot(delx+x0, theta9)
hold off;
legend({'Integral Eq Solution', '1/7th Power Law', '1/9th Power Law'},'Location','northwest');
xlabel('x/L')
ylabel('\theta/L')
print -deps2c Week_2/Exercise4.eps;