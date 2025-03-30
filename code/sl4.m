function srn = sl4(B,Bt,m,n,s)
    k = floor(s/4);
    Om = random_signs(m,k);
    Y = Bt(B(Om));
    Q = orth(Y);
    BQ = B(Q);
    k = s - 3*k;
    Ga = random_signs(n,k);
    Ga = Ga - Q * (Q'*Ga);
    srn = vecnorm(BQ,2,2).^2 + vecnorm(B(Ga),2,2).^2/k;
end