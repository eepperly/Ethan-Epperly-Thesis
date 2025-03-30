function srn = xrownorm_naive(B,Bt,m,n,s)
    k = floor(s/2);
    Om = -3 + 2*randi(2,n,k);
    Y = B(Om);
    srn = zeros(m,1);
    for i = 1:k
        Q = orth(Y(:,[1:(i-1) (i+1):end]));
        C = Bt(Q);
        R = qr(C,"econ");
        z = Q * (Q' * Y(:,i));
        srn = srn + vecnorm(Q*R',2,2).^2 + abs(Y(:,i)-z).^2 ...
            + 2 * real(z .* conj(Y(:,i) - z));
    end
    srn = srn / k;
end