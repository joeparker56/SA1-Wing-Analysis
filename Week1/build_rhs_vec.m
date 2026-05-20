function rhsvec = build_rhs_vec(xs,ys, alpha)
    np = length(xs)-1; 
    rhsvec = zeros(np+1, 1); %column vector b 
    psi_fs = zeros(np,1);

    
    psi_fs(1:np,1) = ys(1:np).*cos(alpha) - xs(1:np).*sin(alpha); %free stream streamfunction
     
   
    rhsvec(1:np-1,1) = psi_fs(1:np-1,1)-psi_fs(2:np,1);
    
    rhsvec(np,1)=0;
    rhsvec(np+1,1)=0;
    %keep b(np) and b(np+1)=0 such that gamma_0=1 and gamma_0 = gamma_np+1
end