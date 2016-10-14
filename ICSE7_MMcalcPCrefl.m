function [MM,cohVals] = ICSE7_MMcalcPCrefl(layer,Lam,AOI,bool_reflect,delta_lam,m_max)

[epsilon,mu,alpha,d,eul] = layer{:};
eul = eul*pi/180; % ZXZ passive Euler rotation angels
AOI = AOI*pi/180; % angle of incidence
MM = zeros(4,4,length(Lam));
A = [1,0,0,1;1,0,0,-1;0,1,1,0;0,1i,-1i,0]; % coherency to Mueller
n0 = 1; % ambient refractive index
Psi0 = [cos(AOI),0,-cos(AOI),0;n0,0,n0,0;0,1,0,1;0,n0*cos(AOI),0,-n0*cos(AOI)]; %ambient medium matrix
Psi0inv = inv(Psi0);
Eul = [cos(eul(3)),sin(eul(3)),0;-sin(eul(3)),cos(eul(3)),0;0,0,1]...
    *[1,0,0;0,cos(eul(2)),sin(eul(2));0,-sin(eul(2)),cos(eul(2))]...
    *[cos(eul(1)),sin(eul(1)),0;-sin(eul(1)),cos(eul(1)),0;0,0,1];
Eul = Eul.';
R = [Eul,zeros(3,3);zeros(3,3),Eul];
kx = n0*sin(AOI); % x component of wave vectors
for n = 1:length(Lam) % inside this loop is the whole Berreman calculation
    M = R*[epsilon(:,:,n),-alpha(:,:,n);transpose(alpha(:,:,n)),mu(:,:,n)]*R.';
    A1 = [M(3,1),M(3,5)+kx;M(6,1),M(6,5)];
    A2 = [M(3,2),-M(3,4);M(6,2)-kx,M(6,4)];
    B1 = [M(5,3)+kx,M(5,6);M(1,3),M(1,6)];
    B2 = [M(4,3),M(4,6);-M(2,3),-M(2,6)+kx];
    C1 = M([3 6],[3 6]);
    D1 = M([5 1],[1 5]);
    D2 = [M(5,2),-M(5,4);M(1,2),-M(1,4)];
    D3 = [M(4,1),M(4,5);-M(2,1),-M(2,5)];
    D4 = [M(4,2),-M(4,4);-M(2,2),M(2,4)];
    temp1 = B1/C1;
    temp2 = B2/C1;
    delta = [-temp1*A1+D1,-temp1*A2+D2;temp2*A1-D3,temp2*A2-D4];
    [Psi1,K] = eig(delta);
    Kdiag = diag(K);
    list1 = find(Kdiag > 0);
    list2 = find(Kdiag < 0);
    if Kdiag(list1(1)) < Kdiag(list1(2))
        list1 = fliplr(list1);
    end
    if Kdiag(list2(1)) > Kdiag(list2(2))
        list2 = fliplr(list2);
    end
    K = diag(exp([Kdiag(list1(1)),Kdiag(list1(2)),Kdiag(list2(1)),Kdiag(list2(2))]...
        *2*pi*d*1i/Lam(n)));
    Psi1 = [Psi1(:,list1(1)),Psi1(:,list1(2)),Psi1(:,list2(1)),Psi1(:,list2(2))];
    P1 = diag([1/K(1,1),1/K(2,2)]);
    P2 = K([3,4],[3,4]);
    invPsi1 = inv(Psi1);
    Int01 = inv(Psi0)*Psi1;
    Int12 = invPsi1*Psi0;
    Int10 = invPsi1*Psi0;
    T01 = inv(Int01([1,2],[1,2]));
    R01 = Int01([3,4],[1,2])*T01;
    T12 = inv(Int12([1,2],[1,2]));
    R12 = Int12([3,4],[1,2])*T12;
    T10 = inv(Int10([3,4],[3,4]));
    R10 = Int10([1,2],[3,4])*T10;
    tau = sum(abs(Kdiag))*d/2;
    
    temp = P2*R12*P1;
    Q = R10*temp;
    D = T10*temp;
    B = T12*P1;
    
    cohVals = ones(1,m_max+1);
    for m = 1:m_max
        cohVals(m+1) = coherence(Lam(n),delta_lam,tau*m);
    end
    if bool_reflect
        Qpowers = zeros(2,2,m_max);
        Qpowers(:,:,1) = eye(2);
        Qpowers(:,:,2) = Q;
        for m = 3:m_max
            Qpowers(:,:,m) = Qpowers(:,:,m-1)*Q;
        end
        temp = zeros(4);
        temp2 = temp;
        temp3 = temp;
        for j = 1:m_max
            temp = temp + cohVals(j+1).*kron(R01,conj(Qpowers(:,:,j)));
            temp2 = temp2 + cohVals(j+1).*kron(Qpowers(:,:,j),conj(R01));
            for k = 1:m_max
                temp3 = temp3 + cohVals(abs(j-k)+1).*kron(Qpowers(:,:,j),conj(Qpowers(:,:,k)));
            end
        end
        coh = kron(R01,conj(R01))...
            + kron(eye(2),conj(D))*temp*kron(eye(2),conj(T01))...
            + kron(D,eye(2))*temp2*kron(T01,eye(2))...
            + kron(D,conj(D))*temp3*kron(T01,conj(T01));
    else
        Qpowers = zeros(2,2,m_max+1);
        Qpowers(:,:,1) = eye(2);
        Qpowers(:,:,2) = Q;
        for m = 3:(m_max+1)
            Qpowers(:,:,m) = Qpowers(:,:,m-1)*Q;
        end
        temp = zeros(4);
        for j = 1:(m_max+1)
            for k = 1:(m_max+1)
                temp = temp + cohVals(abs(j-k)+1).*kron(Qpowers(:,:,j),conj(Qpowers(:,:,k)));
            end
        end
        coh = kron(B,conj(B))*temp*kron(T01,conj(T01));
    end
    MM(:,:,n) = real(A*coh*A')./2;
end
end