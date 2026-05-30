function rhsvec = build_rhs(xs,ys, alpha)
    np = length(xs)-1; 
    rhsvec = zeros(np+1, 1); %column vector b 
    psi_fs = zeros(np,1);

    for i = 1:np
        psi_fs(i,1) = ys(i)*cos(alpha) - xs(i)*sin(alpha); %free stream streamfunction
    end 
    for i = 1:np-1
        rhsvec(i,1) = psi_fs(i,1)-psi_fs(i+1,1);
    end
    rhsvec(np,1)=0;
    rhsvec(np+1,1)=0;
    %keep b(np) and b(np+1)=0 such that gamma_0=1 and gamma_0 = gamma_np+1
end