function srn = sl3(B,Bt,m,n,s)
    k = floor(s/3);
    Om = random_signs(m,k);
    Y = Bt(Om);
    Q = orth(Y);
    BQ = B(Q);
    k = s - 2*k;
    Ga = random_signs(n,k);
    Ga = Ga - Q * (Q'*Ga);
    srn = vecnorm(BQ,2,2).^2 + vecnorm(B(Ga),2,2).^2/k;
end