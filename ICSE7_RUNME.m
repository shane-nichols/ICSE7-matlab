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

% General calculation Parameters, same as in the paper for aBBO, except
% delta_lam which I made smaller here so the plots are more interesting.
d = 13000; % thickness of the crystal, in nanometers
eul = [-8.1,45,0]; % array of passive ZXZ Euler rotation angles, in deg
Lam = 300:750;  % array of measurement wavelengths, in nanometers
AOI = 45; % angle of incidence, in deg
delta_lam = 1.3; % spectral width, in nanometers
bool_reflect = 1; % 1 for reflection, 0 for transmission
fineStep = 0.01; % compute coherent calculation every this many nanometers 
    % for direct integration
n_max = 4; % number of passes through the medium
m_max = 2;  % number of multiple reflections (2 is equivalent to n_max = 4)

% put all the parameters into a cell array
param = {d,eul,Lam,AOI,delta_lam,bool_reflect,fineStep,n_max,m_max};

% function ICSE7_working sets up and performs the calculations 
[MM_dirInt,MM_general,MM_zRecip,MM_multRefl] = ICSE7_working(param);

% % Normalize the calculations 
% % Remove comments here to normalize by M_1,1 (as in the manuscript).
% for n=1:length(Lam)
%     MM_dirInt(:,:,n) = MM_dirInt(:,:,n)./MM_dirInt(1,1,n);
%     MM_general(:,:,n) = MM_general(:,:,n)./MM_general(1,1,n);
%     MM_zRecip(:,:,n) = MM_zRecip(:,:,n)./MM_zRecip(1,1,n);
%     MM_multRefl(:,:,n) = MM_multRefl(:,:,n)./MM_multRefl(1,1,n);
% end

% Below are some examples of using the plotting utility. The program is
% optimized for my 13" MacBook Pro display but it should work on other
% display so long as they are not too difference in size.

% Plot the results with function MMplot. This can be as simple as

% MMplot(Lam,MM_dirInt);

% Below is an more advanced use of the plotter where plots are written
% sequentially. Many things can be configured, such as the size of the
% figure, and the fontsize of all graphics object. A title is added,
% wavelength is converted to eV, a minimum range of the y-axis is set, and
% axes properties are configured. For details, look to the comments in
% MMplot. The third input to MMplot is an optional Matlab LineSpec. So, our
% first call plots MM_dirInt in blue with a linewidth of 2. The handles are
% stored in h.
h = MMplot(Lam,MM_dirInt,'b',...
    'ev',true,...
    'limY',0.05,...
    'size',[1100,700],...
    'fontsize',14,...
    'title','Sequential Plotting example',...
    'lineNV',{'LineWidth',2},...
    'axNV',{...
        'Xgrid','on',...
        'Ygrid','on',...
        'LineWidth',1.5
    }...
);

% Next, plot more lines on the same figure by providing the handles. Most
% options do not need to be handed off again. Exceptions include 'ev' and
% 'limY'. (Note that limY is having no effect here because no MM elements 
% have a range smaller than 0.05, but I show it for example purposes.)
% lineNV is set for each call (obviously as each call plots one line), and
% axNV can be updated for each call (not shown here). Size and fontsize
% should only be set on the first call because the position vectors are not
% recalculted on subsequent calls. This includes resizing by dragging fig
% window! I am working on a callback to permit resizing after the first call.
MMplot(Lam,MM_general,'r',...
    'ev',true,...
    'limY',0.05,...
    'lineNV',{'LineWidth',1.5},...
    'handles',h);

MMplot(Lam,MM_zRecip,'g',...
    'ev',true,...
    'limY',0.05,...
    'lineNV',{'LineWidth',1},...
    'handles',h);

MMplot(Lam,MM_multRefl,'y',...
    'ev',true,...
    'limY',0.05,...
    'lineNV',{'LineWidth',0.5},...
    'handles',h);

% Multiple plots can be written at once. It is better to not set LineSpec
% in this case and let MATLAB pick the colors for you. A legend can be made
% by setting the 'legend' option.

MMplot(Lam,cat(4,MM_dirInt,MM_general,MM_zRecip,MM_multRefl),...
    'ev',true,...
    'limY',0.05,...
    'size',[1100,700],...
    'fontsize',14,...
    'title','Unnormalized Reflection Mueller matrix. Zoom in to see the differences.',...
    'lineNV',{'LineWidth',1.5},...
    'axNV',{...
        'Xgrid','on',...
        'Ygrid','on',...
        'LineWidth',1.5
    },...
    'legend',{'Calculation Type',...
        {'Direct Integration','General Method','z-Reciprocal','Mult. Reflections'}}...
);