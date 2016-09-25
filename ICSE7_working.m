function [MM_dirInt,MM_general,MM_zRecip,MM_multRefl] = ICSE7_working(param)

[d,eul,Lam,AOI,delta_lam,bool_reflect,fineStep,n_max,m_max] = param{:};

% Direct numerical integration. Maybe there is a much better way to do this.
intWidth = 7.2 * delta_lam; % integraion interval
LamLong = (Lam(1)-intWidth/2):fineStep:(Lam(length(Lam))+intWidth/2);

% Optical functions of alpha-barium-borate from www.refractiveindex.info
lam2 = (LamLong/1000).^2;
epsilon(1,1,:) = -lam2*0.0155+0.0184./(lam2 - 0.0179)+2.7405;
epsilon(2,2,:) = epsilon(1,1,:);
epsilon(3,3,:) = -lam2*0.0044+0.0128./(lam2 - 0.0156)+2.373;
mu = repmat(eye(3),[1,1,length(LamLong)]);
alpha = zeros(3,3,length(LamLong));
layer = {epsilon,mu,alpha,d,eul}; % cell array of material parameters

tic
MM_dirInt = ICSE7_MMcalcPCintg(layer,AOI,Lam,bool_reflect,delta_lam,fineStep,intWidth);
disp(['The direct integration finished in ',num2str(toc),' seconds.'])

epsilon = zeros(3,3,length(Lam));
lam2 = (Lam/1000).^2;
epsilon(1,1,:) = -lam2*0.0155+0.0184./(lam2 - 0.0179)+2.7405;
epsilon(2,2,:) = epsilon(1,1,:);
epsilon(3,3,:) = -lam2*0.0044+0.0128./(lam2 - 0.0156)+2.373;
mu = repmat(eye(3),[1,1,length(Lam)]);
alpha = zeros(3,3,length(Lam));
layer = {epsilon,mu,alpha,d,eul}; % cell array of material parameters

tic
MM_general = ICSE7_MMcalcPCperm(layer,Lam,AOI,bool_reflect,delta_lam,n_max);
disp(['The general calculation finished in ',num2str(toc),' seconds.'])

tic
MM_zRecip = ICSE7_MMcalcPCcomb(layer,Lam,AOI,bool_reflect,delta_lam,n_max);
disp(['The z-reciprocal calculation finished in ',num2str(toc),' seconds.'])

tic
MM_multRefl = ICSE7_MMcalcPCrefl(layer,Lam,AOI,bool_reflect,delta_lam,m_max);
disp(['The multiple reflections calculation finished in ',num2str(toc),' seconds.'])

end

