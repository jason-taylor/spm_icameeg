function ICA = spm_icameeg_scorr(S)
% Correlate ICA weight topographies with artefact topographies.
%  FORMAT: ICA = spm_icameeg_scorr(S)
%  INPUT: Struct 'S' with fields:
%   S.D          - MEEG object or filename of MEEG object
%   S.Dcont      - Filename of continuous data (if S.D is epoched)
%   S.ICA        - ICA struct or filename (output of spm_eeglab_runica)
%   S.artchans   - artefact chan names, cell array (default: {'VEOG'})
%   S.arttypes   - types of artefacts to find (spm_eeg_artefact_?) 
%                    (default: {'eyeblink'})
%   S.artmax     - threshold to reject artefact trials (default: 1000) uV
%   S.thresh     - threshold for zscore of spatial correlations (def: 2)
%   S.rmfiles    - 1/0 remove intermediate files? (def: 1)
%  OUTPUT: 
%   ICA          - ICA struct with artefact.scorr field with subfields: 
%    (arttype).rmat    - summary matrix with columns: 
%                         component index, r, z-score of r
%    (arttype).topo    - vector summary of artefact topography
%    (arttype).suspect - indices of suprathreshold components
%
%  NOTE: Currently only works for 'eyeblink'.
%
%  spm_icameeg and spm_eeglab tools
%  by Jason Taylor (09/Mar/2018) jason.taylor@manchester.ac.uk
%   + 18/May/2021 jt: save fig as well as png, fstem ICA.fname not D.fname
%   + 20/Jun/2025 jt: added artefact rejection of extreme-valued trials

%--------------------------------------------------

%% Check inputs:
if ~isstruct(S.ICA)
    load(S.ICA);
else
    ICA = S.ICA;
    S.ICA = ICA.fname; % save RAM
end
try artchans = S.artchans; catch, artchans = {'VEOG'};     end
try arttypes = S.arttypes; catch, arttypes = {'eyeblink'}; end
try thresh = S.thresh;     catch, thresh = 2;              end
try fnamecont = S.Dcont;   catch, fnamecont = '';          end
try rmfiles = S.rmfiles;   catch, rmfiles = 1;             end
try artmax = S.artmax;     catch, artmax = 1000;           end


%% Load SPM-format data file (on which ICA was run):

Dica = spm_eeg_load(S.D);

% Load continuous data (if D is epoched):
if ~isempty(fnamecont)
    D = spm_eeg_load(fnamecont);
else
    D = Dica;
end

fprintf('\n\n');
fprintf('++ %s\n',datestr(now));
fprintf('++ RUNNING spm_icameeg_scorr ON %s\n',D.fname);
arttypetxt = sprintf('%s ',arttypes{:});
fprintf('++ USING: artypes: %s\n',arttypetxt);
artchantxt = sprintf('%s ',artchans{:});
fprintf('++ USING: artchans: %s\n',artchantxt);
fprintf('++ USING: thresh =  %g\n',thresh);

% Get inverse of ICA weights for topographies:
iweights = pinv(ICA.weights);

files2delete = {};

% Prepare figure:
spm_figure('Clear','Graphics');
fig = spm_figure('GetWin','Graphics');
if ~isnumeric(fig)
    fig=fig.Number;
end

for i=1:length(arttypes)
    artchan = artchans{i};
    arttype = arttypes{i};
    suspect = [];
     
    %% Filter (for artefact finding):
    S = [];
    S.D     = D.fname;
    S.band  = 'bandpass';
    S.freq  = [1 20];
    S.type  = 'butterworth';
    S.order = 5;
    S.dir   = 'twopass';
    
    Do = spm_eeg_filter(S);
    files2delete = [files2delete;Do.fname;Do.fnamedat];
    
    %% Find artefacts using spm_eeg_artefact:
    S = [];
    S.D                          = Do.fname;
    S.mode                       = 'mark';
    S.badchanthresh              = 0.2;
    S.methods.channels           = artchan;
    S.methods.fun                = arttype;
    S.methods.settings.threshold = 4;
    S.methods.settings.excwin    = 200; % was 0, but that now fails
    S.append                     = true;
    S.prefix                     = 'art_tmp_';
    
    Do = spm_eeg_artefact(S);
    files2delete = [files2delete;Do.fname;Do.fnamedat];
    
    %% Epoch (artefacts):
    S = [];
    S.D                       = Do.fname;
    S.timewin                 = [-100 300];
    S.trialdef.conditionlabel = arttype;
    S.trialdef.eventtype      = sprintf('artefact_%s',arttype);
    S.trialdef.eventvalue     = artchan;
    S.trialdef.trlshift       = 0; % was -200
    S.bc                      = 1;
    S.prefix                  = 'e_';
    S.eventpadding            = 0;
    
    Do = spm_eeg_epochs(S);
    files2delete = [files2delete;Do.fname;Do.fnamedat];

    %% Reject Artefacts - remove events with extreme values
    S = [];
    S.D = Do.fname;
    S.mode = 'reject';
    S.badchanthresh = 0.2;
    S.prefix = 'a';
    S.append = true;
    S.methods(1).channels = artchan;
    S.methods(1).fun = 'threshchan';
    S.methods(1).settings.threshold = artmax; % user set or 1000uV default
    S.methods(1).settings.excwin = 1000;

    Do = spm_eeg_artefact(S);
    files2delete = [files2delete;Do.fname;Do.fnamedat];
    
    %% Average (artefacts):
    S = [];
    S.D           = Do.fname;
    S.robust      = false;
    S.circularise = false;
    S.prefix      = 'm';
    
    Do = spm_eeg_average(S);
    files2delete = [files2delete;Do.fname;Do.fnamedat];
    
    %% Get artefact distribution over channels:
    chanlabs = ICA.chans;
    chaninds = indchannel(Do,ICA.chans);
    c2d = coor2D(Do,chaninds);

    d = selectdata(Do,chanlabs,[.050 .250],[]);
    dm = mean(d,2);
    
    % Plot topo:
    spm_eeg_plotScalpData(dm,c2d,chanlabs);
    clim = get(gca,'clim');
    colormap jet
    cm = colormap;
    oldax = gca; oldfig = gcf;
    spm_figure('Clear','Graphics');
    spm_figure('Focus',fig);
    sp = subplot(3,1,1);
    axis ij equal off;
    copyobj(get(oldax,'children'),sp);
    title(arttype);
    set(gca,'clim',clim);
    cm(1,:) = [1 1 1];
    colormap(cm);
    set(findobj(gca,'Type','Line'),'visible','off');
    close(oldfig)

    %% Compute spatial correlations:
    rs = zeros(ICA.ncomp,1);
    for j=1:length(rs)
        try
            rs(j) = abs(corr(dm,iweights(:,j)));
        catch
            tmp = corrcoef(dm,iweights(:,j));
            rs(j) = abs(tmp(2));
        end
    end
    
    % Compute z-score of correlations (for thresholding)
    zrs = (rs-mean(rs))/std(rs);

    % Construct summary matrix, sort by z score (descending)
    rmat = [(1:length(rs))' rs zrs];
    rmat = sortrows(rmat,2);

    %% Plot z-scores:
    subplot(3,1,3); 
    bar(zrs); title(arttype);
    xlabel('IC'); ylabel('z-score of correlation');
    set(gca,'xtick',1:ICA.ncomp);
    set(gca,'fontsize',8);
    axis tight

    % Report suprathreshold correlations:
    indsupthresh = find(rmat(:,3)>thresh);
    fprintf('++ Found %d suprathreshold spatial correlations with %s\n',length(indsupthresh),artchan);
    if any(indsupthresh)
        suspect(1,:) = rmat(indsupthresh,1)';
        for k=indsupthresh'
            fprintf('IC%d:\tr=%0.4f\tz=%0.4f\n',rmat(k,:));
        end
    else
        fprintf('++ Highest 3 correlations:\n');
        fprintf('IC%d:\tr=%0.4f\tz=%0.4f\n',rmat(end-2,:));
        fprintf('IC%d:\tr=%0.4f\tz=%0.4f\n',rmat(end-1,:));
        fprintf('IC%d:\tr=%0.4f\tz=%0.4f\n',rmat(end,:));
    end
    
    % Plot suspect component topography(ies):
    if any(suspect)
        S=[];
        S.D              = Dica.fname;
        S.ICA            = ICA;
        S.components     = suspect;
        S.figsubplot     = [fig 3 length(suspect) length(suspect)+1];
        S.savefig        = 0;
        spm_icameeg_topo(S);
    end
    
    % Store results in ICA:
    ICA.artefact.scorr.(arttype).summary = rmat;
    ICA.artefact.scorr.(arttype).suspect = suspect;
    ICA.artefact.scorr.(arttype).thresh  = thresh;
    ICA.artefact.scorr.(arttype).topo    = dm;

    %% Save figure:
    [~,fstem] = fileparts(ICA.fname);
    figfname = sprintf('summary_scorr_%s_%s.png',arttype,fstem);
    print(fig,'-dpng',figfname);
    fprintf('++ Saved figure to image: %s\n',figfname);
    figfname = sprintf('summary_scorr_%s_%s.fig',arttype,fstem);
    saveas(fig,figfname,'fig');
    fprintf('++ Saved figure to fig: %s\n',figfname);
    
end

%% Save, delete files:

% Save ICA:
save(ICA.fname,'ICA');
fprintf('++ Saved output ICA struct to %s\n',ICA.fname);

% Delete intermediate files:
if rmfiles
    fprintf('++ Deleting intermediate files\n');
    for i=1:length(files2delete)
        delete(files2delete{i})
    end
end

return