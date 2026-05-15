function rhsvec = build_rhs(xs,ys, alpha)
    np = length(xs)-1;
    rhsvec = zeros(np+1, 1);
    psi_fs = zeros(np,1);

    for i = 1:np
        psi_fs(i,1) = ys(i)*cos(alpha) - xs(i)*sin(alpha);
    end 
    for i = 1:np-1
        rhsvec(i,1) = psi_fs(i,1)-psi(i+1,1);
    end
end