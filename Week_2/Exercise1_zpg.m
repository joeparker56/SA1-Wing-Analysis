clear;
close all;

np = 101;
x = linspace(0,1,np);
u_ = 1*ones(1,np); %dimensionless velocity 

int = 0; % value of the integral 
theta_2 = zeros(1,np);  %theta squared at each position 
re_l = 1000;
for i=2:np
    int = int + ueintbit(x(i-1), u_(i-1), x(i), u_(i));
    theta_2(i) = (0.45/re_l) * u_(i)^(-6) * int;
end
analy = 0.664/re_l^(0.5) * x.^0.5;
theta = sqrt(theta_2);


figure;
plot(x, theta);
hold on;
plot(x, analy);
hold off;
xlabel('x/L');
ylabel('\theta/L');
legend('Computed \theta', 'Blausius Solution');
grid on;
print -deps2c Week_2/Exercise1.eps;
