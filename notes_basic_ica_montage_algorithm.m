% Basic algorithm for applying ICA to SPM12-format EEG data via montage.

% Load SPM data:
D = spm_eeg_load; % select continuous data file

% Extract data matrix (or could substitute for d in runica call)
% NOTE: if epoched data, use 'reshape' to coerce to 2D
chtype = 'EEG';
d = selectdata(D,chanlabels(D,indchantype(D,chtype)),[],[]);

% Run ICA
npc = 32; % number of PCA dimensions
[weights,sphere,compvars,bias,signs,lrates,activations] = runica(d,'extended',1,'pca',npc);

% Decide which components to remove
% NOTE: Could look for correlations between activations and EOG/ECG.
% NOTE: Could look for artefact topographies in channel weight topos.
comp2remove = [1 7]; % for example

% Create SPM montage matrix w/ inverse ICA weights, removing artefact ICs
iweights = pinv(weights);
finalics = setdiff(1:size(weights,1),comp2rem);
tramat = iweights(:,finalics) * ICA.weights(finalics,:);

montage = [];
montage.tra = tramat;
montage.labelorg = chanlabels(D,indchantype(D,chtype));
montage.labelnew = chanlabels(D,indchantype(D,chtype));

save('montage_ICA_cleaned_EEG.mat','montage');

% Now apply this montage to your data.
% Or I suppose you could run this separately for each channel type, then
% combine the montages into one montage file and apply that. 
