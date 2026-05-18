function psixy = psipv(xc, yc, Gamma, x, y)
    r_2 = (x-xc).^2 + (y-yc).^2;
    psixy = - (Gamma / (4 * pi)) * log(r_2);
end

