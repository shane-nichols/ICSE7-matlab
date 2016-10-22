function handles = MMplot(Lam,MMdata,varargin)
% Mueller matrix 2D plotting utility
% Shane Nichols, New York University

% Makes a 4 x 4 array of 2-D line plots with full control over line and
% axes properties.
% Outputs: [1 x 16] array of axis handles
%
% Required positional inputs:
%   Lam: [1 x n] array of wavelengths (X-axis)
%   MMdata: [4 x 4 x n x ...] Mueller matrix array
% Optional positional inputs:
%   LineSpec: string containing a valid lineSpec. Type "doc LineSpec" in
%       command window for more info. Default is "-", a solid line.
% Optional Name-Value pairs inputs:
%   ev: bool. converts X axis to eV. e.g., 'ev',true
%   handles: [1 x 16] array of plot handles. New handles are created if not given.
%   limY: scalar numeric. limits how small the range of the y-axes can be.
%   fontsize: sets font-size. Default is 12 pts. Changing the fontsize
%       of existing plots is not recommended. (Set on first call).
%   lineNV: a 1D cell array containing Name-Value pair arguments valid for
%       Chart Line Properties.
%   axNV: a 1D cell array containing Name-Value pairs arguments valid for
%       Axes Properties.
%   size: Size of the figure in pixels given as a two element vector [X Y].
%       A warning is issued if the requested size is larger than the screen
%       size minus the height of the OSX status bar (on my machine). 
%       Default size is [1000 700].
%   title: string containing a title to place at the top of the figure.
%   legend: two-element cell array. First element is a string to use for
%       title of the legend. Second element is either a numeric array 
%       containing values to use for labels of each plot, or a cell array
%       of strings to use as labels. Only set legend on last call, or just
%       write all plots at once (better).

p = inputParser;
% input validation functions
valFun1 = @(x) ischar(x) && ...
    all(~strcmpi(x,{'ev','handles','lineNV','limY','fontsize','axNV','size','title','legend'}));
valFun2 = @(x) isscalar(x)&&isnumeric(x);
% setup input scheme
addRequired(p,'Lam',@isnumeric);
addRequired(p,'MMdata',@isnumeric);
addOptional(p,'LineSpec','-',valFun1)
addParameter(p,'ev',false,@islogical)
addParameter(p,'handles',gobjects(1,16), @(x) all(ishandle(x)))
addParameter(p,'limY',0,valFun2)
addParameter(p,'fontsize',12,valFun2)
addParameter(p,'axNV',{},@iscell)
addParameter(p,'lineNV',{},@iscell)
addParameter(p,'size',[1000 700],@(x) length(x) == 2 && isnumeric(x))
addParameter(p,'title','',@ischar)
addParameter(p,'legend',{},@iscell)
parse(p,Lam,MMdata,varargin{:}) %parse inputs

% create new figure if no valid handles were given
handles = p.Results.handles;
if any(strcmpi('handles',p.UsingDefaults))
    % Determine how large to make the figure window, according to the screensize.
    scrsz = get(0,'screensize');
    figPos = [1 5 p.Results.size];
    if figPos(3) > scrsz(3)
        figPos(3) = scrsz(3);
        warning(['Figure horizontal dimension set to the maximum value of ',...
            num2str(figPos(3)),' pixels.'])
    end
    if figPos(4) > (scrsz(4) - 99)   % 99 pixels is the height of the OSX status bar on my machine
        figPos(4) = (scrsz(4) - 99);
        warning(['Figure vertical dimension set to the maximum value of ',...
            num2str(figPos(4)),' pixels.'])
    end
    h_fig = figure('position',figPos,'units','pixels'); %create figure
    xLabel = uicontrol('style','text','BackgroundColor','w',...
        'units','pixels','FontSize',p.Results.fontsize); % create x-label
    if p.Results.ev == true
        set(xLabel,'String','Energy (eV)');
    else
        set(xLabel,'String','Wavelength (nm)');
    end
    xLabel_sz = get(xLabel,'extent');
    set(xLabel,'Position',[(figPos(3) - xLabel_sz(3) )./2, 0, xLabel_sz(3), xLabel_sz(4)]);
    
    if ~isempty(p.Results.title) % create title if given
        figTitle = uicontrol('style','text','BackgroundColor','w',...
            'units','pixels','FontSize',p.Results.fontsize);
        set(figTitle,'String',p.Results.title)
        figTitle_sz = get(figTitle,'extent');
        set(figTitle,'Position',[( figPos(3) - figTitle_sz(3) )./2,...
            ( figPos(4) - figTitle_sz(4) ), figTitle_sz(3), figTitle_sz(4)]);
    end
    % determine the horizontal extent of y-axis marker labels
    dummy = uicontrol('style','text','fontsize',p.Results.fontsize,'units','pixels');
    set(dummy,'String','-0.000');
    yAxSz = get(dummy,'extent');
    delete(dummy)
    
    plotSzX = figPos(3)/4 - yAxSz(3) - yAxSz(3)./5; % X size of plot area in pixels
    plotSzY = ( figPos(4) - 4*yAxSz(4) )/4 - 6; % Y size of plot area in pixels
    for i=1:4
        for j=1:4
            plotPos = [ ( (plotSzX + yAxSz(3) + 3)*(j-1) + yAxSz(3) +5)./figPos(3) ,...
                ((plotSzY + yAxSz(4)./2)*(4-i)+yAxSz(4)*2 + 3)./figPos(4),...
                plotSzX./figPos(3), plotSzY./figPos(4)];
            hand = subplot('Position',plotPos);
            hold(hand,'on')
            box(hand,'on')
            if i ~= 4
                set(hand,'XTickLabel',[]) % keep X lables only for bottom row
            end
            handles(j+4*(i-1)) = hand;
        end
    end
else
    h_fig = get(handles(1),'parent');
    figPos = get(h_fig,'Position');
end

%plot data and set Line properties.
if p.Results.ev == true; Lam = 1239.8./Lam; end
for j = 1:4
    for k = 1:4
        plot(handles(k+4*(j-1)),Lam,squeeze(MMdata(j,k,:,:)),...
            p.Results.LineSpec,p.Results.lineNV{:})
    end
end
% set Axes properties
axis(handles,'tight'); % first, axes are set to tight
if ~isempty(p.Results.axNV)
    for j=1:16; set(handles(j),p.Results.axNV{:}); end
end
if p.Results.limY ~= 0 % modify axes bounds if limY is set
    lim = p.Results.limY;
    for j=1:16
        Ylim = get(handles(j),'YLim');
        if (Ylim(2) - Ylim(1)) < lim
            avg = (Ylim(2) + Ylim(1))./2;
            Ylim(2) = avg + lim/2;
            Ylim(1) = avg - lim/2;
            set(handles(j),'Ylim',Ylim);
        end
    end
end
% Adjust plot limits so that lines do not overlap axis borders. 
% *** If you like to use Markers, then perhaps change 'lineWidth' to 'MarkerSize'
lineHandle = get(handles(1),'children');
lineWidth = zeros(size(lineHandle));
for j = 1:length(lineHandle)
    lineWidth(j) = get(lineHandle(j),'lineWidth');
end
lineWidth = max(lineWidth);
plotPos = get(handles(1),'Position');
for j=1:16
    xlim = get(handles(j),'xLim');
    ylim = get(handles(j),'yLim');
    xStep = (xlim(2) - xlim(1))/plotPos(3)/figPos(3)*lineWidth/2;
    yStep = (ylim(2) - ylim(1))/plotPos(4)/figPos(3)*lineWidth;
    set(handles(j),'XLim',[xlim(1)-xStep,xlim(2)+xStep]);
    set(handles(j),'YLim',[ylim(1)-yStep,ylim(2)+yStep]);
end
% set font size of all graphics objects if fontsize was passed
if ~any(strcmpi('fontsize',p.UsingDefaults))
    set(get(gcf,'children'),'FontSize',p.Results.fontsize);
end
% optionally create legend (this will increase the width of the figure!)
if ~any(strcmpi('legend',p.UsingDefaults))
    Labels = p.Results.legend{2};
    if isnumeric(Labels); Labels = strread(num2str(Labels),'%s'); end
    for i=1:16
        set(handles(i),'units','pixels');
        pos(:,i) = get(handles(i),'Position');
    end
    lgd = legend(handles(4),Labels,'location','northeastoutside');
    set(lgd,'units','pixels','fontsize',p.Results.fontsize);
    title(lgd,p.Results.legend{1},'FontSize',p.Results.fontsize);
    lgd_pos = get(lgd,'Position');
    h_fig.Position = h_fig.Position + [0 0 lgd_pos(3) 0];
    for i=1:16
        set(handles(i),'Position',pos(:,i));
    end
end

end