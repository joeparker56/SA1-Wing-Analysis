function [theta, Hes, int, ils, itr, its,x] = calc_bl(dudx, re_l)
    global re_l dudx u0

    np = 101;
    %for turbulent seperation at x/L=1, dudx = -0.383
    x = linspace(0,1,np); %length varies from 0->1
    
        
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
        if Hes(i)>1.58 && ils ~= 0 &&itr==0
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
        Hes(its+1:np) = Hes(i);
        theta(its+1:np) = theta(i).*(u_(i)./u_(i+1:np)).^(H+2);
    end

end
