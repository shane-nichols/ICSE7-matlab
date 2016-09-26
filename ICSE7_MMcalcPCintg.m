function MM_DirInt = ICSE7_MMcalcPCintg(layer,AOI,Lam,bool_reflect,delta_lam,fineStep,intWidth)

LamLong = (Lam(1)-intWidth/2):fineStep:(Lam(length(Lam))+intWidth/2);
[epsilon,mu,alpha,d,eul] = layer{:};
eul = eul*pi/180; % ZXZ passive Euler rotation angels
AOI = AOI*pi/180; % angle of incidence
MM = zeros(4,4,length(LamLong));
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
for n = 1:length(LamLong) % inside this loop is the whole Berreman calculation
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
    [Psi,K] = eig(delta); % eig is faster than expm (although this could be made analytic)
    P = Psi0inv*Psi*diag(exp(diag(K)*2*pi*d*1i/LamLong(n)))/Psi*Psi0;
    if bool_reflect
        J =  P([3 4],[1 2])/(P([1 2],[1 2]));
        MM(:,:,n) =  real(A*kron(J,conj(J))*A'./2);
    else
        J = inv(P([1 2],[1 2]));
        MM(:,:,n) =  real(A*kron(J,conj(J))*A'./2);
    end
end

gaussian = exp(-(-(intWidth/2):fineStep:(intWidth/2)).^2./(2*delta_lam^2));
sumGauss = sum(gaussian);
lengthGauss = length(gaussian);
if mod(lengthGauss,2)
    halfWidth = lengthGauss/2+1/2;
else
    halfWidth = lengthGauss/2+1;
end
idx = 0;
MM_DirInt = zeros(4,4,length(Lam)); % do summations at values in Lam.
idx2 = round(fracIndex(LamLong,Lam));
for n=Lam
    idx = idx+1;
    temp=zeros(4);
    for k=1:lengthGauss
        temp = temp + MM(:,:,idx2(idx)-halfWidth+k).*gaussian(k);
    end
    MM_DirInt(:,:,idx)=temp./sumGauss; % normalization
end
end