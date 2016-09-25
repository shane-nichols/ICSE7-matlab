function MM = ICSE7_MMcalcPCcomb(layer,Lam,AOI,bool_reflect,delta_lam,n_max)
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
for lam = 1:length(Lam) % inside this loop is the whole Berreman calculation
    M = R*[epsilon(:,:,lam),-alpha(:,:,lam);transpose(alpha(:,:,lam)),mu(:,:,lam)]*R.';
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
    q = [Kdiag(list1(1)),Kdiag(list1(2)),Kdiag(list2(1)),Kdiag(list2(2))]; % sorted eigenvalues
    K = diag(exp(q.*2*pi*d*1i./Lam(lam)));
    Psi1 = [Psi1(:,list1(1)),Psi1(:,list1(2)),Psi1(:,list2(1)),Psi1(:,list2(2))];
    P1 = diag([1/K(1,1),1/K(2,2)]);
    P2 = K([3,4],[3,4]);
    invPsi1 = inv(Psi1);
    Int01 = Psi0inv*Psi1;
    Int12 = invPsi1*Psi0;
    Int10 = invPsi1*Psi0;
    T01 = inv(Int01([1,2],[1,2]));
    R01 = Int01([3,4],[1,2])*T01;
    T12 = inv(Int12([1,2],[1,2]));
    R12 = Int12([3,4],[1,2])*T12;
    T10 = inv(Int10([3,4],[3,4]));
    R10 = Int10([1,2],[3,4])*T10;
    % recursive calculation of wave combinations
    if bool_reflect
        n_max = fix(n_max/2)*2; % make sure n_max is even for reflection
    else
        n_max = fix((n_max+1)/2)*2-1; % make sure n_max is odd for transmission
    end
    U = zeros(2,2,(n_max+1)*(n_max+2)/2-1);
    U(:,:,1) = [P1(1,1)*T01(1,1) , P1(1,1)*T01(1,2) ; 0 , 0]; % Eq. 17a
    U(:,:,2) = [0 ,0 ; P1(2,2)*T01(2,1) , P1(2,2)*T01(2,2)];  % Eq. 17b
    for n=2:n_max % this loop does Eqs. 18 and 19
        idx = (n+1)*n/2;
        if mod(n,2)
            U(:,:,idx) = [P1(1,1)*R10(1,1)*U(1,:,idx-n) ; 0,0];
            for k=1:(n-1)
                U(:,:,idx+k) = [P1(1,1)*R10(1,:)*U(:,:,idx-n+k) ;...
                    P1(2,2)*R10(2,:)*U(:,:,idx-n+k-1)];
            end
            U(:,:,idx+n) = [0,0 ; P1(2,2)*R10(2,2)*U(2,:,idx-1)];
        else
            U(:,:,idx) = [P2(1,1)*R12(1,1)*U(1,:,idx-n) ; 0,0];
            for k=1:(n-1)
                U(:,:,idx+k) = [P2(1,1)*R12(1,:)*U(:,:,idx-n+k) ;...
                    P2(2,2)*R12(2,:)*U(:,:,idx-n+k-1)];
            end
            U(:,:,idx+n) = [0,0 ; P2(2,2)*R12(2,2)*U(2,:,idx-1)];
        end
    end
    if bool_reflect == false % this block will make the array in Eq. 20a
        n2 = fix((n_max+1)/2);
        n2 = n2.*(n2+1);
        J = zeros(2,2,n2);
        times = zeros(1,n2);
        idx2 = 1;
        for n = 1:2:n_max
            idx = (n+1)*n/2;
            for k=0:n
                J(:,:,idx2) = U(:,:,idx+k);
                times(idx2) = d*(k*q(2)+(n-k)*q(1));
                idx2 = idx2 + 1;
            end
        end
    else % this block will make array in Eq. 20b
        n2 = fix((n+2)/2).^2;
        J = zeros(2,2,n2);
        times = zeros(1,n2);
        idx2 = 1;
        for n = 2:2:n_max
            idx = (n+1)*n/2;
            for k=0:n
                J(:,:,idx2) = U(:,:,idx+k);
                times(idx2) = d*(k*q(2)+(n-k)*q(1));
                idx2 = idx2 + 1;
            end
        end
        J(:,:,idx2) = inv(T10)*R01;
        times(:,:,idx2) = 0;
    end
    coh = zeros(4);
    for i=1:n2  % this block does the double sum in Eq. 22
        for j=1:n2
            tau = times(i)-times(j);
            coh = coh + coherence(Lam(lam),delta_lam,tau)*kron(J(:,:,i),conj(J(:,:,j)));
        end
    end
    if bool_reflect % propagate the field into the ambient
        coh = kron(T10,conj(T10))*coh;
    else
        coh = kron(T12,conj(T12))*coh;
    end
    MM(:,:,lam) = real(A*coh*A')./2; %convert to Mueller
end