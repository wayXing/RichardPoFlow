function [phi,Uu,P] = DEIM(U)

    [n,m]=size(U);

    [rho,phi_1] = max(abs(U(:,1)));
    Uu=U(:,1);
    P = zeros(n,1); 
    P(phi_1,1) = 1;
    phi=phi_1;
    
    for i=2:m
        c=(P'*Uu)\(P'*U(:,i));
        r=U(:,i)-Uu*c;
        [rho,phi_i] = max(abs(r));
        
        Uu = [Uu,U(:,i)];
        P_i = zeros(n,1); P_i(phi_i,1)=1;
        P = [P,P_i];
        phi = [phi; phi_i];
    end
end