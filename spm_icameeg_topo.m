function ftopo = spm_icameeg_topo(S)
% Plot ICA channel-weight topographies (from spm_eelab_runica).
%  FORMAT: ftopo = spm_icameeg_topo(S)
%  INPUT: Struct 'S' with fields:
%   S.D          - MEEG object or filename of MEEG object
%   S.ICA        - ICA struct or filename (output of spm_eeglab_runica)
%   S.components - Indices of components to plot (default: all)
%   S.figsubplot - [fig sprows spcols startind] to move plot to
%   S.savefig    - save figure to file? (1(or prefix)/0 def: 1)
%  OUTPUT: 
%   ftopo        - Figure handle of topography plot
%
%  spm_icameeg and spm_eeglab tools
%  by Jason Taylor (09/Mar/2018) jason.taylor@manchester.ac.uk
%   + 18/May/2021 jt: save plot to image and figure; chan type
%   + 05/Aug/2021 jt: added prefix option for topo fig/image filename
%
%--------------------------------------------------

%% Check inputs:
if ~isstruct(S.ICA)
    load(S.ICA);
else
    ICA = S.ICA;
    S.ICA = ICA.fname; % save RAM
end
try components = S.components; catch, components = 1:ICA.ncomp; end
try figsubplot = S.figsubplot; catch, figsubplot = [];          end
try savefig    = S.savefig;    catch, savefig    = 1;           end
    
%% Load SPM-format data file, get channel info:
D = spm_eeg_load(S.D);
cl = ICA.chans;
c2d  = D.coor2D(D.indchannel(cl));
if length(cl)==1
    cl = input('Cell array of chanlabels for topo?: ');
end

fprintf('\n\n');
fprintf('++ %s\n',datestr(now));
fprintf('++ RUNNING spm_icameeg_topo ON %s\n',ICA.fname);

% Inverse weight matrix:
iweights = pinv(ICA.weights);
    
%% Setup figure/subplots
if isempty(figsubplot)
    spm_figure('Clear','Graphics');
    ftopo = spm_figure('GetWin','Graphics');
    if ~isnumeric(ftopo)
        ftopo = ftopo.Number;
    end
    if length(components)<5
        spdims = [length(components) 1];
    elseif mod(sqrt(length(components)),1)==0
        spdims = [1 1]*sqrt(length(components));
    else
        r = floor(sqrt(length(components))); % jt 3/3/2018 - fixed error
        c = ceil(length(components)/floor(sqrt(length(components)))); %jt fixed
        spdims = [r c];
    end
    spstart = 1;
else
    ftopo   = figsubplot(1);
    spdims  = figsubplot(2:3);
    spstart = figsubplot(4);
end
    
% hack: prevent backwards text!
try opengl('software');
catch, feature('UseGenericOpenGL',1);
end

%% Plot topographies

% options for scalp maps:
in=[];
in.noButtons = 1;
in.cbar = 0;
in.plotpos = 0;
in.type = char(chantype(D,indchannel(D,cl{1})));

for i=1:length(components)
    comp = components(i);
    
    % Plot topo:
    spm_eeg_plotScalpData(iweights(:,comp),c2d,cl,in);
    oldax = gca; oldfig = gcf;
    colormap jet; % jt 3/3/2018 - def colormap is different in new matlab
    cm = colormap;
    clim = get(oldax,'clim');
    
    % Copy to summary figure: 
    figure(ftopo);
    sp = subplot(spdims(1),spdims(2), i+spstart-1);
    axis ij equal off;
    copyobj(get(oldax,'children'),sp);
    title(sprintf('IC%02d',comp));
    close(oldfig)
    set(gca,'clim',clim);

end

% Format figure:
cm(1,:) = [1 1 1];
colormap(cm)
set(gcf,'color','w')

%% Save figure as image:
if savefig
    [~,fstem] = fileparts(ICA.fname);
    if isstr(savefig)
        figfstem = sprintf('%s_%s',savefig,fstem);
    else
        figfstem = sprintf('topos_%s',fstem);
    end
    
    pngfname = [figfstem '.png'];
    print(gcf,'-dpng',pngfname);
    fprintf('++ Figure saved to image: %s\n',pngfname);

    figfname = [figfstem '.fig'];
    saveas(gcf,figfname,'fig');
    fprintf('++ Figure saved to fig: %s\n',figfname);

end

return
