clear;
close all;


global re_l dudx u0

re_l = 1e7;
dudx =  -0.5;
u0 = 1;

x0 = 0.01;
thick0(1) = 0.023*x0*(re_l*x0)^(-1/6);
thick0(2) = 1.83*thick0(1);

[delx, thickhist] =ode45(@thickdash, [0, 0.99], thick0);

x=x0+delx;

all_He = thickhist(:,2) ./ thickhist(:,1);

sep_location=[];

for i = 1:length(all_He)
    if all_He(i) < 1.46 && isempty(sep_location)
        sep_location = i;
    end
end


disp(['Re_l ', num2str(re_l)]);
disp(['dudx ', num2str(dudx)]);
if ~isempty(sep_location)
    disp(['Separation at x/L = ', num2str(x(sep_location), '%.4f')]);
else
    disp('No separation before x/L = 1');
end



figure;
plot(x, thickhist(:,1), 'DisplayName', '\theta/L')
hold on;
plot(x,thickhist(:,2), 'DisplayName', '\delta_e/L' )
xlabel('x/L'); ylabel('\theta/L,  \delta_E/L');
legend;
print -deps2c Week_2/Exercise5.eps;


figure;
plot(x, all_He);
xlabel('x/L'); ylabel('He');
yline(1.46, '--','Separation threshold');
