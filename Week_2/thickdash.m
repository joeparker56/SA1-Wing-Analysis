
function dthickdx = thickdash(xmx0, thick)
    global Re dudx u0

    dthickdx = zeros(2,1);
    He = thick(2)/thick(1);
    if He>=1.46
        H =(11*He + 15)/(48*He-59);
    else
        H = 2.803;        
    end 
    ue = u0 + dudx*xmx0;
    re_thet = Re * ue * thick(1);
    cf = 0.091416*((H-1)*re_thet).^(-0.232) * exp(-1.26*H);
    cdiss = 0.010026*((H-1)*re_thet).^(-1/6);
    dthickdx(1,1) = cf/2 - ((H+2)/ue)*dudx*thick(1);
    dthickdx(2,1) = cdiss - (3/ue)*dudx*thick(2);
end