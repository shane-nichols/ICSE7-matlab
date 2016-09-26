% This script computes a Mueller matrix array for alpha-Barium borate using
% four differnet approaches, as explained in the paper:
%   MM_dirInt:   by direct integration (Eq. 1)
%   MM_general:  by wave permutations (Eqs. 14)
%   MM_zRecip:   by wave combinations (Eqs. 22)
%   MM_multRefl: by partially coherent multiple reflections (Eqs. 23)

% The calculation times will be printed to the command window. Because the
% general and the z-reciprocal calculations require many nested for-loops,
% they are not very fast in MatLab. We expect they would be much faster in
% a lower level language like C++.

% General calculation Parameters
d = 13000; % thickness of the crystal, in nanometers
eul = [-8.1,45,0]; % array of passive ZXZ Euler rotation angles, in deg
Lam = 300:750;  % array of measurement wavelengths, in nanometers
AOI = 45; % angle of incidence, in deg
delta_lam = 1.3; % spectral width, in nanometers
bool_reflect = 1; % 1 for reflection, 0 for transmission
fineStep = 0.01; % compute coherent calculation every this many nanometers for direct integration
n_max = 4; % number of passes through the medium
m_max = 2;  % number of multiple reflections to include

param = {d,eul,Lam,AOI,delta_lam,bool_reflect,fineStep,n_max,m_max};

[MM_dirInt,MM_general,MM_zRecip,MM_multRefl] = ICSE7_working(param);

% % normalize the calculations. 
% Remove comments below to normalize by M_1,1 (as in the manuscript).
% for n=1:length(Lam)
%     MM_dirInt(:,:,n) = MM_dirInt(:,:,n)./MM_dirInt(1,1,n);
%     MM_general(:,:,n) = MM_general(:,:,n)./MM_general(1,1,n);
%     MM_zRecip(:,:,n) = MM_zRecip(:,:,n)./MM_zRecip(1,1,n);
%     MM_multRefl(:,:,n) = MM_multRefl(:,:,n)./MM_multRefl(1,1,n);
% end

% plot the results with linewidths getting smaller so overlap can be seen.
h = MMplot(Lam,MM_dirInt,'-b',...
    'ev',true,...
    'title','Unnormalized Reflection Mueller matrix',...
    'lineNV',{'LineWidth',2});

MMplot(Lam,MM_general,'-r',...
    'ev',true,...
    'title','Unnormalized Reflection Mueller matrix',...
    'lineNV',{'LineWidth',1.5},...
    'handles',h);

MMplot(Lam,MM_zRecip,'-g',...
    'ev',true,...
    'title','Unnormalized Reflection Mueller matrix',...
    'lineNV',{'LineWidth',1},...
    'handles',h);

MMplot(Lam,MM_multRefl,'-y',...
    'ev',true,...
    'title','Unnormalized Reflection Mueller matrix',...
    'lineNV',{'LineWidth',0.5},...
    'handles',h);

% Discrepency between the direct integration and the other methods is
% proportional to the derivative of the optical functions.
