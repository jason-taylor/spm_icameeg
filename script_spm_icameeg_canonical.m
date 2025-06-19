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
% RUN each cell one at a time (click to highlight cell, then Ctrl+Enter)
%
% by Jason Taylor (18/Jun/2019) jason.taylor@manchester.ac.uk
%  + jt (19/Jun/2025) updated interpolation to work with SPM25


%% Load SPM12 M/EEG format data

% Epoched: Will run faster, waste fewer components on between-trial noise.
% Continuous: Slower, but sometimes easier to interpret component 
%  time-courses (because you see everything, not just what happened during
%  a trial). Can waste components modelling between-trial noise, which
%  might be a problem if you have large gaps between trials.

D = spm_eeg_load;
fname = D.fname;
[~,fstem] = fileparts(fname);
fprintf('\nFilename is %s\n\n',fname);


%% Plot data to view interactively (with EOG, ECG data linked)

S=[];
S.D              = D.fname;
S.chantype       = 'EEG'; % 'MEGMAG','MEGPLANAR'
S.title          = 'Original EEG Data';
S.spacing        = 100;
S.child.chantype = {'EOG'};
S.child.title    = 'EOG';
S.child.spacing  = 300;

spm_eeglab_eegplot(S);


%% Create bad-channel interpolation montage

% Requires hacked function ft_channelrepair_jt.m to be copied here:
%  [spm path]/external/fieldtrip/
% Uncomment and run code below... (should now work with SPM25 too)

% if strmatch('SPM25',spm('version'))
%     src = which('COPYME_ft_channelrepair_jt_SPM25.m');
% else
%     src = which('COPYME_ft_channelrepair_jt.m');
% end
% spmdir = fileparts(which('spm'));
% dst = fullfile(spmdir,'external','fieldtrip','ft_channelrepair_jt.m');
% copyfile(src,dst);

S=[];
S.D            = D.fname;
S.montagefname = sprintf('montage_badchans_%s.mat',fstem);
S.fixbads      = 0; % 1 if you want to interpolate before ICA
                    % 0 if you want to omit bad chans from ICA and
                    %   interpolate them later (during clean-ICA montage)

[D,montagefname_bc,montage_bc] = spm_uomeeg_channelrepair(S);


%% Run ICA (takes a few minutes)

% This will return a struct 'ICA' with fields containing output from
% EEGLAB's 'runica' and some other information. Functions that plot IC
% topographies or timecourses or look for artefact components will require
% both this ICA struct and the SPM12 M/EEG format data. ICA will be updated
% by some of these functions.

S=[];
S.D        = D.fname;
S.chantype = 'EEG'; % only use EEG (EOG will still contain artefact)
S.ncomp    = 32; % number of components (PCA before ICA)
S.skipbad  = 1; % 1=don't use bad channels (will interpolate later)
S.prefix   = sprintf('ICA%d_',S.ncomp); % <- eg. 'ICA32_'
S.seed     = 0; % to be able to replicate (seed for random # generator)

ICA = spm_eeglab_runica(S);


%% Correlate IC time-courses with artefact (EOG) channels

% If any component has captured the blink signal, we would expect it to
% have a high temporal correlation with VEOG. An eye-movement component
% might correlate with both HEOG and VEOG, depending on how participants
% move their eyes. I flag components with a particulalry high correlation
% (based on z-score r-mean(r)/std(r) over component correlations).

S=[];
S.D        = D.fname;
S.ICA      = ICA.fname;
S.artchans = {'VEOG','HEOG'};
%S.artchans = {'EEG062','EEG061'};
S.artfilts = {[1 20],[1 20]};
S.thresh   = 2; % z-score of correlations; default=2;

ICA = spm_icameeg_tcorr(S);


%% Correlate IC weights with artefact topographies

% If a component captures an artefact, then its channel weights should have
% a similar topographic distribution to the artefact itself. SPM has a
% routine that identifies and marks blinks in continuous data. Here I've
% done that, then epoched and averaged around the blink events to compute
% an average blink topography. I then correlate that topography to each
% component's vector of channel weights and flag any with particularly high
% correlations (based on z-score r-mean(r)/std(r) over component 
% correlations).

S=[];
S.D        = D.fname; % if epoched, use S.Dcont below...
S.Dcont    = S.D(3:end);
%S.Dcont    = S.D(2:end); % if fname epoched; assumes D.fname is e*.mat
S.ICA      = ICA.fname;
S.artchans = {'VEOG'};
%S.artchans = {'EEG062'};
S.arttypes = {'eyeblink'};
S.thresh   = 2; % z-score of spatial correlations; default=2
S.rmfiles  = 1; % remove intermediate files (blink epochs)

ICA = spm_icameeg_scorr(S);


%% Plot ICA channel weight topographies

% Note: If you identify any suspect ICs, you can store them in the ICA
% struct using, e.g. (increment '1' as needed),
%   ICA.artefact.manual(1).suspect = 13;
%   ICA.artefact.manual(1).reason  = 'muscle artefact topo?';
%   save(ICA.fname,'ICA');

S=[];
S.D          = D.fname;
S.ICA        = ICA.fname;
S.components = [1:ICA.ncomp]; % <- eg. [1 2 3 ... 32]

spm_icameeg_topo(S);


%% Plot ICA 'activations' (with EOG, ECG data linked)

% Note: If you identify any suspect ICs, you can store them in the ICA
% struct using, e.g. (increment '1' as needed),
%   ICA.artefact.manual(1).suspect = 13;
%   ICA.artefact.manual(1).reason  = 'high-frequency EMG artefact?';
%   save(ICA.fname,'ICA');

S=[];
S.D              = D.fname;
S.ICA            = ICA.fname;
S.compinds       = 1:ICA.ncomp;
S.title          = 'IC Activations';
S.spacing        = 5;
S.child.chantype = {'EOG'};
S.child.title    = 'EOG';
S.child.spacing  = 300;

spm_eeglab_eegplot(S);


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

S=[];
S.D         = D.fname;
S.ICA       = ICA.fname;
S.dotopo    = 1; % plot topographies of suspect components?
S.doeegplot = 1; % plot time-courses of suspect components?

[suspects,report] = spm_icameeg_suspects(S);


%% Create montage to remove artefact ICs (and interpolate bad chans)

% By now you should have determined which components you want to remove
% from the data. Enter their numbers in S.remove below, or use comp2remove
% below to select components labelled 'suspect' automatically.

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


%% Plot 'cleaned' data (with EOG, ECG data linked)

S=[];
S.D              = D.fname; % <- ICA-cleaned data
S.D2             = fname;   % <- original 'dirty' data (plotted in red)
S.chanlabs       = ICA.chans;
S.title          = 'Clean(blue) and Original(red) Data';
S.spacing        = 100;
S.child.chantype = {'EOG'};
S.child.title    = 'EOG';
S.child.spacing  = 300;

spm_eeglab_eegplot(S);


%% Plot topographies of removed components and report artefact info

% Save topography image for later reference. Report correlations etc. Can
% be useful for compiling information for group studies (e.g., number of
% components removed, average VEOG correlation of removed components, etc.)

S=[];
S.D         = D.fname;
S.ICA       = ICA.fname;
S.dotopo    = 1; % plot topographies of removed components?
S.savetopo  = 'topo_removed'; % prefix for .fig and .png
S.doeegplot = 0; % plot time-courses of removed components?

[removed,report] = spm_icameeg_removed(S);


%% Continue preprocessing...

% You can now carry on with your usual preprocessing -- but don't reject 
% based artefact channels if you have removed the corresponding component 
% (e.g., VEOG if you have removed blinks)!
