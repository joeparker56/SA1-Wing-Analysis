clear;
close all;

global re_l dudx u0

np = 101;

x = linspace(0,1,np); %length varies from 0->1
dudx = -0.5; %dimensionless gradient
re_l = 1e5;

theta = zeros(1,np);  %theta squared at each position
Hes = zeros(1,np);
Hes(1) = 1.57258;


laminar = true;
integ = 0;
int = 0; %location of laminar->turbulent transition
ils = 0; %location of laminar seperation
itr = 0; %location of turbulent reattachment 
its = 0; %location of turbulent seperation
i=1;
u_ = linspace(1, 1+dudx, np);
while laminar && i<np
    i=i+1;
    integ = integ + ueintbit(x(i-1), u_(i-1), x(i), u_(i));
    theta(i) = sqrt((0.45/re_l) * u_(i)^(-6) * integ);
    m = -re_l * theta(i)^2  * dudx;
    H = thwaites_lookup(m);
    Hes(i) = laminar_He(H);
    re_thet = re_l * u_(i) * theta(i);
    if log(re_thet) >= 18.4*Hes(i) - 21.74
        laminar = false;    %transition
        int = i;

    elseif m>= 0.09
        laminar = false;
        ils = i;
        Hes(i) = 1.51509;
    end            
end
delta_e = Hes(i)*theta(i);

while its == 0 && i<np %loop for turbulent boundary layer
    i=i+1;
    thick0 = [theta(i-1), delta_e];
    u0 = u_(i-1); %velocity at start of panle
    [delx, thickhist] = ode45(@thickdash, [0, x(i)-x(i-1)], thick0);
    theta(i) = thickhist(end, 1);
    delta_e = thickhist(end, 2);
    Hes(i) = delta_e/theta(i);
    if Hes(i)>1.58 && ils ~= 0
        itr = i; %turbulent reattachment
    end
    if Hes(i) <1.46
        its = i; %turbulent separation
    end
end

if Hes(i)>=1.46
        H =(11*Hes(i) + 15)/(48*Hes(i)-59);
else
    H = 2.803;        
end 
if its ~= 0
    Hes(i+1:np) = Hes(i);
    theta(i+1:np) = theta(i).*(u_(i)./u_(i+1:np)).^(H+2);
end


disp(['re_l = ' num2str(re_l) ', dudx = ' num2str(dudx)])
if int ~= 0
    disp(['Natural transition at x = ' num2str(x(int)) ', Rethet/1000 = ' num2str(re_l * u_(int) * theta(int) / 1000)])
end
if ils ~= 0
    disp(['Separation at x = ' num2str(x(ils)) ', Rethet/1000 = ' num2str(re_l * u_(ils) * theta(ils) / 1000)])
end
if itr ~= 0
    disp(['Turbulent Reattachment at x = ' num2str(x(itr)) ', Rethet/1000 = ' num2str(re_l * u_(itr) * theta(itr) / 1000)])
end
if its ~= 0  
    disp(['Turbulent Seperation at x = ' num2str(x(its)) ', Rethet/1000 = ' num2str(re_l * u_(its) * theta(its) / 1000)])
end
if int==0 && ils ==0
    disp('No transition or separation')
end
disp(' ')


figure;
plot(x, theta);
xlabel('x/L');
ylabel('\theta/L');
if int ~= 0
    xline(x(int), '-', 'Laminar Transition')
end
if ils ~= 0
    xline(x(ils), '-', 'Laminar Seperation')
end
if itr ~= 0
    xline(x(itr), '-', 'Turbulent Reattachment')
end
if its ~= 0  
    xline(x(its), '-', 'Turbulent Seperation')
end

figure;
plot(x, Hes);
xlabel('x/L');
ylabel('He')
if int ~= 0
    xline(x(int), '-', 'Laminar Transition')
end
if ils ~= 0
    xline(x(ils), '-', 'Laminar Seperation')
end
if itr ~= 0
    xline(x(itr), '-', 'Turbulent Reattachment')
end
if its ~= 0  
    xline(x(its), '-', 'Turbulent Seperation')
end
