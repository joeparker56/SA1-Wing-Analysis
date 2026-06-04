function [int ils itr its delstar theta] = bl_solv(x,cp)
    global Re dudx u0
    np = length(x);
    delstar = zeros(1,np);
    ue = sqrt(1-cp);
    ue_ext = [1,ue];
    x_ext = [0,x];
    dudxs = (ue_ext(2:np)-ue_ext(1:np-1))./(x_ext(2:np)-x_ext(1:np-1));

    theta = zeros(1,np);
    delstar = zeros(1,np);
    Hes = zeros(1,np);
    laminar = true;

    integ = 0;

    int = 0;
    ils = 0;
    itr = 0;
    its = 0;

    dudx = dudxs(1);
    integ = ueintbit(0,0, x(1), ue(1));
    theta(1) = sqrt((0.45/Re) * ue(1)^(-6) * integ);
    m = - Re * theta(1)^2  * dudx;
    H = thwaites_lookup(m);
    delstar(1) = H*theta(1);
    Hes(1) = laminar_He(H);
    i=1;

    while laminar && i<np
        i=i+1;
        dudx = dudxs(i);
        integ = integ + ueintbit(x(i-1), ue(i-1), x(i), ue(i));
        theta(i) = sqrt((0.45/Re) * ue(i)^(-6) * integ);
        m = - Re * theta(i)^2 * dudx;
        H = thwaites_lookup(m);
        delstar(i) = H * theta(i);
        Hes(i) = laminar_He(H);
        re_thet = Re*ue(i) * theta(i);
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

    while its ==0 && i<np
        i=i+1;
        thick0 = [theta(i-1), delta_e];
        u0=ue(i-1);
        dudx = dudxs(i-1);
        [delx, thickhist] = ode45(@thickdash, [0,x(i)-x(i-1)], thick0);
        theta(i) = thickhist(end, 1);
        delta_e = thickhist(end, 2);
        Hes(i) = delta_e/theta(i);
        if Hes(i)>=1.46
            H =(11*Hes(i) + 15)/(48*Hes(i)-59);
        else
            H = 2.803;        
        end 
        delstar(i) = H*theta(i);        
        if Hes(i)>1.58 && ils ~= 0 &&itr==0
            itr = i; %turbulent reattachment
        end
        if Hes(i) <1.46
            its = i; %turbulent separation
        end
    end
    if its~=0
        for i = its+1:np
            Hes(i)=Hes(its);
            if Hes(i)>=1.46
                H =(11*Hes(i) + 15)/(48*Hes(i)-59);
            else
                H = 2.803;        
            end 
            theta(i) = theta(its)*(ue(its)/ue(i))^(H+2);
            delstar(i) = H*theta(i);
        end
    end
end


