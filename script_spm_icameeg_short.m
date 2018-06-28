% Short, automated ICA script using 'spm_icameeg', 'spm_eeglab', and 
% 'spm_uomeeg' tools.
%  by Jason Taylor (09/Mar/2018) jason.taylor@manchester.ac.uk

%% Data:

cd C:\work\data\myproject\subject01
fname  = 'ffMdspmeeg_subject01.mat'; % epoched or pre-epoching file
%cfname = 'ffMdspmeeg_subject01.mat'; % if fname epoched (for scorr)

%% PATHS:

addpath C:\work\software\spm12_git\
addpath C:\work\software\eeglab_v14_git\
addpath C:\work\software\spm_icameeg
addpath C:\work\software\spm_eeglab
addpath C:\work\software\spm_uomeeg

spm eeg;
eeglab; 

%% Run ICA:

S=[];
S.D              = fname;
S.ncomp          = 32;
S.chantype       = 'EEG';
S.skipbad        = 1;  % <- requires badchan interp montage (later)
S.seed           = 0;
S.prefix         = 'ICA_';

ICA = spm_eeglab_runica(S);

%% Correlate time-courses with VEOG, HEOG:

S=[];
S.D              = fname;
S.ICA            = ICA.fname;
S.artchans       = {'VEOG','HEOG'};
S.artfilts       = {[1 20],[1 20]};
S.thresh         = 2;

ICA = spm_icameeg_tcorr(S);

%% Correlate channel weights with blink topography:

S=[];
S.D              = fname; % if epoched, use S.Dcont below...
%S.Dcont          = cfname;
S.ICA            = ICA.fname;
S.artchans       = {'VEOG'};
S.arttypes       = {'eyeblink'};
S.thresh         = 2; % z-score of spatial correlations; default=2
S.rmfiles        = 1; % remove intermediate files (blink epochs)

ICA = spm_icameeg_scorr(S);

%% Report suspect components:

S=[];
S.D              = fname;
S.ICA            = ICA.fname;
S.dotopo         = 1;
S.doeegplot      = 0;

[suspects,report] = spm_icameeg_suspects(S);

%% Create an interpolation montage for any bad channels:

S=[];
S.D              = fname;
S.montagefname   = 'montage_badchans.mat';

[D,montagefname_bc,montage_bc] = spm_uomeeg_channelrepair(S);

%% Project artefact components out of data:

S=[];
S.D              = fname;
S.ICA            = ICA.fname;
S.remove         = suspects;
S.montagefname   = 'montage_ICA_clean.mat';
S.bcmontagefname = montagefname_bc;
S.apply          = 1;

[D,ICA,montagefname,montage] = spm_icameeg_montage(S);
