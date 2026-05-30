function lhsmat = build_lhs_fully_vec(xs, ys)
    np = length(xs)-1;
    
    %XA are the starting points of our vortex-sheet-panels
    %X_eval are the points at which we will evaluate the streamfunction
    %However, we will evaluate the streamfunction at every starting point
    %of our sheet panels, so we create two matrices with pairs of both at each
    %index
    [XA, X_eval] = meshgrid(xs(1:np), xs(1:np)');
    [YA, Y_eval] = meshgrid(ys(1:np), ys(1:np)');
    
    %These are the end points of each vortex-sheet-panel. Each index of the
    %endpoints corresponds to the start point for the same panel, so we can
    %use element-wise calculations.
    [XB,~] = meshgrid(xs(2:np+1), xs(1:np)');
    [YB, ~] = meshgrid(ys(2:np+1), ys(1:np)');
    
    %Create a matrix of the influence coefficients. 
    %The row of the matrix corresponds to the point on which its acting.
    %The column of the matrix corresponds to the sheet panel from which it
    %comes from. 

    [infa, infb] = panelinf_vec(XA, YA, XB, YB, X_eval, Y_eval);

    %Here will fill in the psip matrix based on the influence coefficients
    %at each point.
    psip = zeros(np, np+1);
    psip(:,1) = infa(:,1);
    psip(:,2:np) = infa(:,2:np) + infb(:,1:np-1);
    psip(:, np+1) = infb(:,np);
    
    %We set up the A matrix by finding the difference between the i+1th and
    %ith item in the psip matrix
    lhsmat = zeros(np+1, np+1);
    lhsmat(1:np-1, :) = psip(2:np,:)-psip(1:np-1,:);

    %set up the equation that gamma_np+1 = gamma_1
    lhsmat(np, 1) = 1;
    lhsmat(np, np+1) = -1;
    
    %set up the equation that gamma_1 = 0
    lhsmat(np+1,1)=-1;
    lhsmat(np+1,2)=+1;
    lhsmat(np+1,3)=-0.5;
    lhsmat(np+1,np-1)=-1;
    lhsmat(np+1,np)=0.5;
end
