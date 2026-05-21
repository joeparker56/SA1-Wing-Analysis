clear;
close all;

np = 101;
x = linspace(0,1,np);
u_ = 1*ones(1,np); %dimensionless velocity 

int = 0; % value of the integral 
theta_2 = zeros(1,np-1);  %theta squared at each position 
re_l = 1000;
for i=1:np-1
    int = int + ueintbit(x(i), u_(i), x(i+1), u_(i+1));
    
    theta_2(i) = (0.45/re_l) * u_(i)^(-6) * int;
end 

theta = sqrt(theta_2);
plot(x(1:np-1), theta);
