% Script to run ICA on SPM12 M/EEG data using 'spm_icameeg', 'spm_eeglab', 
% and 'spm_uomeeg' tools.
% 
% ICA is run using EEGLAB's 'runica'.
% Components are plotted (time-courses using 'eegplot' and channel-weight
% topographies using spm_eeg_plotScalpData) and correlated with artefact
% time-courses or topographies (using custom functions).
% Artefact components are projected out of the data using spm_eeg_montage.
%
% Custom functions call EEGLAB and SPM12 functions (and some FieldTrip
% functions that are bundled with SPM12).
%
% * SET 'fname' and 'datadir' (and toolbox paths) at the top.
% * SET processing flags to execute/skip each step.
% * RUN each cell one at a time (click to highlight cell, then Ctrl+Enter)
%    or run it all start to finish if you're brave!
%
%  by Jason Taylor (09/Mar/2018) jason.taylor@manchester.ac.uk

%% Parameters & Paths -- EDIT THESE AS NECESSARY!

% Data (EDIT THESE!):
datadir = 'path/to/your/data';
fname   = 'my_spmeeg_data.mat'; % if epoched, use continuous below:
%cfname  = 'my_continuous_spmeeg_data.mat'; % ony needed if fname epoched

% Toolbox paths (EDIT THESE!):
spm_path    = 'C:\work\software\spm12_git';
eeglab_path = 'C:\work\software\eeglab_git';
spm_icameeg_path = 'C:\work\software\spm_icameeg';
spm_eeglab_path  = 'C:\work\software\spm_eeglab';
spm_uomeeg_path  = 'C:\work\software\spm_uomeeg';

% processing flags 1=run/0=don't (*=required, ~=recommended, o=optional):
do_plot_cleandata   = 1; % o view original data in 'eegplot'
do_findbadchans     = 0; % o find bad channels (low corr with neighbours)
do_badchan_montage  = 0; % o create badchan montage (will apply later)
do_runica           = 1; % * RUN ICA
do_tcorr_eog        = 1; % ~ temporal correlations IC-EOG
do_tcorr_eeg        = 1; % ~ temporal correlations IC-EEG
do_scorr_eyeblink   = 1; % o spatial correlation w/ blink topo (auto detected)
do_plot_icatopos    = 1; % ~ view IC channel weight topographies
do_plot_activations = 1; % o view IC activations in 'eegplot'
do_report_suspects  = 1; % ~ report automatically labelled suspect components
do_ica_montage      = 1; % * CREATE MONTAGE (& apply it) to remove components
do_plot_cleandirty  = 1; % o view clean & original ('dirty') in 'eegplot'

% Debug option (on error, will stop and give access to function stack):
dbstop if error 
%dbclear if error % removes the above


%% Set paths, start SPM and EEGLAB (if you haven't already)

% SPM:
if isempty(which('spm'))
    addpath(spm_path);
end
spmmenufig = spm_figure('FindWin','Menu');
if isempty(spmmenufig) || isempty(get(spmmenufig,'Name'))
    spm eeg
end

% EEGLAB:
if isempty(which('eeglab'))
    addpath(eeglab_path);
end
if isempty(which('eegplot'))
    eeglab rebuild
end

% Add toolbox paths:
addpath(spm_icameeg_path);
addpath(spm_eeglab_path);
addpath(spm_uomeeg_path);


%% Load the data to start with

cd(datadir);
D = spm_eeg_load(fname);
[~,fstem] = fileparts(fname);


%% (optional) Plot data (with EOG, ECG data linked)

if do_plot_cleandata
    
    S=[];
    S.D              = D.fname;
    S.chantype       = 'EEG'; % 'MEGMAG','MEGPLANAR'
    S.title          = 'Original EEG Data';
    S.spacing        = 100;
    S.child.chantype = {'EOG'};
    S.child.title    = 'EOG';
    S.child.spacing  = 300;
    
    spm_eeglab_eegplot(S);
    
end

%% (optional) Find bad channels based on low correlations with neighbours

if do_findbadchans
    
    S=[];
    S.D            = D.fname;
    S.thresh       = -3; % default=-3 (z-score)
    S.nNN          = 3; % default=3 (number of nearest neighbours)
    
    [D,bads] = spm_uomeeg_findbadchans(S);

    % If it goes wrong, unset using:
    %  D = badchannels(D,bads,0);
    %  D.save;
    
end

%% (optional) Create bad-channel interpolation montage

if do_badchan_montage

    S=[];
    S.D            = D.fname;
    S.montagefname = sprintf('montage_badchans_%s.mat',fstem);
    S.fixbads      = 0; % 1 if you want to interpolate before ICA...
    %S.newprefix    = 'Mbc_'; % only used if fixbads=1
    
    [D,montagefname_bc,montage_bc] = spm_uomeeg_channelrepair(S);

end

%% Run ICA (takes a few minutes)

% This will return a struct 'ICA' with fields containing output from
% EEGLAB's 'runica' and some other information. Functions that plot IC
% topographies or timecourses or look for artefact components will require
% both this ICA struct and the SPM12 M/EEG format data. ICA will be updated
% by some of these functions.

if do_runica 
    
    S=[];
    S.D        = D.fname;
    %S.timewin  = D.time([1 end]); % default (same if left out)
    %S.timewin  = [10 600]; % specify a time-window
    %S.timewin  = {[10 400] [405 600]}; % ... or several
    S.chantype = 'EEG'; % only use EEG (EOG will still contain artefact)
    S.ncomp    = 32; % number of components (PCA before ICA)
    %S.chantype = 'MEGPLANAR';
    %S.ncomp     = 64;
    S.skipbad  = 1; % 1=don't use bad channels (will interpolate later)
    S.prefix   = sprintf('ICA%d_',S.ncomp); % <- eg. 'ICA32_'
    S.seed     = 0; % to be able to replicate (seed for random # generator)
    
    ICA = spm_eeglab_runica(S);

end

%% (recommended) Correlate IC time-courses with artefact (EOG) channels

% If any component has captured the blink signal, we would expect it to
% have a high temporal correlation with VEOG. An eye-movement component
% might correlate with both HEOG and VEOG, depending on how participants
% move their eyes. I flag components with a particulalry high correlation
% (based on z-score r-mean(r)/std(r) over component correlations).

if do_tcorr_eog
    
    S=[];
    S.D        = D.fname;
    S.ICA      = ICA.fname;
    S.artchans = {'VEOG','HEOG'};
    S.artfilts = {[1 20],[1 20]};
    %S.artchans = {'EEG062','EEG061','EEG063'};
    %S.artfilts = {[1 20],[1 20],[1 48]};
    S.thresh   = 2; % z-score of correlations; default=2;
    
    ICA = spm_icameeg_tcorr(S);

end

%% (recommended) Correlate IC time-courses with all EEG channel time-courses

% If any component has captured the noise from a single bad channel, then
% its temporal correlation with that channel will be high whereas its
% correlations with other channels will be low. Here we look for components
% whose ratio of maximum:next-highest channel correlations is very high.

if do_tcorr_eeg

    S=[];
    S.D      = D.fname;
    S.ICA    = ICA.fname;
    S.thresh = 10; % max:next-highest corr ratio; default=10 
    
    ICA = spm_icameeg_chancorr(S);

end

%% (optional) Correlate IC weights with artefact topographies

% If a component captures an artefact, then its channel weights should have
% a similar topographic distribution to the artefact itself. SPM has a
% routine that identifies and marks blinks in continuous data. Here I've
% done that, then epoched and averaged around the blink events to compute
% an average blink topography. I then correlate that topography to each
% component's vector of channel weights and flag any with particularly high
% correlations (based on z-score r-mean(r)/std(r) over component 
% correlations).

if do_scorr_eyeblink
    
    S=[];
    S.D        = D.fname; % if epoched, use S.Dcont below...
    %S.Dcont    = cfname; % if you've specified it, and if fname is epoched.
    %S.Dcont    = S.D(2:end); % if fname epoched; assumes D.fname is e*.mat
    S.ICA      = ICA.fname;
    S.artchans = {'VEOG'};
    %S.artchans = {'EEG062'};
    S.arttypes = {'eyeblink'};
    S.thresh   = 2; % z-score of spatial correlations; default=2
    S.rmfiles  = 1; % remove intermediate files (blink epochs)
    
    ICA = spm_icameeg_scorr(S);
    
end

%% (optional) Plot ICA channel weight topographies

% Note: If you identify any suspect ICs, you can store them in the ICA
% struct using, e.g. (increment '1' as needed),
%   ICA.artefact.manual(1).suspect = 13;
%   ICA.artefact.manual(1).reason  = 'muscle artefact topo?';
%   save(ICA.fname,'ICA');

if do_plot_icatopos
    
    S=[];
    S.D          = D.fname;
    S.ICA        = ICA.fname;
    S.components = [1:ICA.ncomp]; % <- eg. [1 2 3 ... 32]
    
    spm_icameeg_topo(S);

end

%% (optional) Plot ICA 'activations' (with EOG, ECG data linked)

% Note: If you identify any suspect ICs, you can store them in the ICA
% struct using, e.g. (increment '1' as needed),
%   ICA.artefact.manual(1).suspect = 13;
%   ICA.artefact.manual(1).reason  = 'high-frequency EMG artefact?';
%   save(ICA.fname,'ICA');

if do_plot_activations
    
    S=[];
    S.D              = D.fname;
    S.ICA            = ICA.fname;
    S.compinds       = 1:ICA.ncomp;
    S.title          = 'IC Activations';
    S.spacing        = 5;
    S.child.chantype = {'EOG'};
    S.child.title    = 'EOG';
    %S.child.chanlabs = {'EEG061','EEG062','EEG063'};
    %S.child.title    = 'EOG and ECG';
    S.child.spacing  = 300;
    
    spm_eeglab_eegplot(S);
    
end

%% Manually add any more 'suspect' components 

% This bit of code will ask you to manually indicate any suspect components
% and give a reason for each.

i=0;
addmore = input('Do you want to add suspect components manually? (y/n): ','s');
while strcmpi('y',addmore)
    i=i+1;
    ICA.artefact.manual(i).suspect = input('Suspect component: ');
    ICA.artefact.manual(i).reason  = input('Reason suspected: ','s');
    addmore = input('Add more? (y/n): ','s');
end

save(ICA.fname,'ICA');


%% Summary - report all suspect components

if do_report_suspects
    
    S=[];
    S.D         = D.fname;
    S.ICA       = ICA.fname;
    S.dotopo    = 1; % plot topographies of suspect components?
    S.doeegplot = 1; % plot time-courses of suspect components?
    
    [suspects,report] = spm_icameeg_suspects(S);

end

%% Create montage to remove artefact ICs (and interpolate bad chans)

% By now you should have determined which components you want to remove
% from the data. Enter their numbers in S.remove below, or use comp2remove
% below to select components labelled 'suspect' automatically.

if do_ica_montage

    % Automatic decision:
    %comp2remove = suspects; % accept decisions from above

    % Interactive decision (with prompt):
    fprintf('Suspects are: [ %s]\n',sprintf('%d ',suspects));
    comp2remove = input('Which components to remove? (include [ ] for multiple): ');
    fprintf('REMOVING [ %s]\n',sprintf('%d ',comp2remove));
    
    S=[];
    S.D              = D.fname;
    S.ICA            = ICA.fname;
    S.remove         = comp2remove; % <- from above
    S.bcmontagefname = montagefname_bc; % if you haven't interpolated bad chans!
    S.apply          = 1; % apply it to clean the data 
    
    [D,ICA,montagefname,montage] = spm_icameeg_montage(S);

end

%% (optional) Plot 'cleaned' data (with EOG, ECG data linked)

if do_plot_cleandirty
    
    S=[];
    S.D              = D.fname; % <- ICA-cleaned data
    S.D2             = fname;   % <- original 'dirty' data (plotted in red)
    S.chanlabs       = ICA.chans;
    S.title          = 'Clean(blue) and Original(red) Data';
    S.spacing        = 100;
    S.child.chantype = {'EOG'};
    S.child.title    = 'EOG';
    %S.child.chanlabs = {'EEG061','EEG062','EEG063'};
    %S.child.title    = 'EOG and ECG';
    S.child.spacing  = 300;
    
    spm_eeglab_eegplot(S);
    
end

S = [];


%% Continue preprocessing...

% You can now carry on with your usual preprocessing -- but don't reject 
% based artefact channels if you have removed the corresponding component 
% (e.g., VEOG if you have removed blinks)!
