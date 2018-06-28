function [D,ICA,montagefname,montage] = spm_icameeg_montage(S)
% Create a montage to remove ICA components (from spm_eeglab_runica).
%  FORMAT: [ICA,montage,montagefname] = spm_icameeg_montage(S)
%  INPUT: Struct 'S' with fields:
%   S.D              - MEEG object or filename of MEEG object
%   S.ICA            - ICA struct or filename (output of spm_eeglab_runica)
%   S.remove         - Indices of components to remove (def: prompt)
%   S.montagefname   - Filename for output montage
%   S.bcmontagefname - Filename of bad-chan interpolation montage (optional)
%   S.apply          - 1/0 Apply resulting montage? (def: 0 (no))
%   S.newprefix      - (if apply=1) Output prefix for 'clean' data file
%  OUTPUT: 
%   D                - Cleaned data (if apply=1)
%   ICA              - With 'removed' and 'montagefname' subfields
%   montage          - Montage created
%   montagefname     - Filename of montage
%
%  spm_icameeg and spm_eeglab tools
%  by Jason Taylor (09/Mar/2018) jason.taylor@manchester.ac.uk

%--------------------------------------------------

%% Check inputs:
if ~isstruct(S.ICA)
    load(S.ICA);
else
    ICA = S.ICA;
    S.ICA = ICA.fname; % save RAM
end
try 
    comp2rem = S.remove;
catch
    comp2rem = input('components to remove?: ');
end
try bcmontagefname = S.bcmontagefname; catch, bcmontagefname = []; end
try apply = S.apply;                   catch, apply = 0;           end
try newprefix = S.newprefix;           catch, newprefix = 'MICA_'; end

%% Load SPM-format data file:
D = spm_eeg_load(S.D);
[~,fstem] = fileparts(ICA.fname);
try montagefname = S.montagefname; catch, montagefname = []; end

fprintf('\n\n');
fprintf('++ %s\n',datestr(now));
fprintf('++ RUNNING spm_icameeg_montage ON %s\n',D.fname);
if ~isempty(bcmontagefname)
    fprintf('++ USING: bad-channel montage =  %s\n',bcmontagefname);
end
if apply
    fprintf('++ USING: apply montage with prefix: %s\n',newprefix);
end

%% Compute ICA projection montage:

% Inverse weight matrix:
iweights = pinv(ICA.weights);
finalics = setdiff(1:ICA.ncomp,comp2rem);
tramat = iweights(:,finalics) * ICA.weights(finalics,:);

% Put into full montage (all channels, not just EEG):
tra = eye(size(D,1));
tra(D.indchannel(ICA.chans),D.indchannel(ICA.chans)) = tramat;

% Interpolate bad channels?:
msg=''; bctxt='';
if any(bcmontagefname)
    load(bcmontagefname);
    trabc = montage.tra;
    tra = trabc * tra;
    msg = 'and interpolating bad channels';
    bctxt = 'bc_';
end

% Create SPM-format montage:
clear montage
montage.tra = tra;
montage.labelorg = D.chanlabels;
montage.labelnew = D.chanlabels;
if isempty(montagefname)
    remtxt = sprintf('%d_',comp2rem);
    montagefname = sprintf('montage_ICA%d_rem%s%s%s.mat',ICA.ncomp,remtxt,bctxt,fstem);
end
save(montagefname,'montage');
fprintf('++ Created montage: %s\n',montagefname);

%% Plot scaled image:
%f=figure('color','w');
spm_figure('Clear','Graphics');
f = spm_figure('GetWin','Graphics');
imagesc(tra);
colormap jet
axis image
title(['Montage after removing artefact components' msg]);
set(gca,'xtick',1:size(tra,1));
set(gca,'ytick',1:size(tra,2));
set(gca,'fontsize',8);

% Save figure as image:
figfname = sprintf('imagesc_%s.png',montagefname(1:end-4));
print(f,'-dpng',figfname);
fprintf('++ Saved figure to %s\n',figfname);

%% Apply montage?:

if apply 
    
    fprintf('++ Applying montage to data\n');
    
    S=[];
    S.D             = D.fname;
    S.mode          = 'write';
    S.prefix        = newprefix;
    S.montage       = montagefname;
    S.keepothers    = 0; % Must be 0 in case badchan interp is used!
    S.keepsensors   = 1;
    S.updatehistory = 1;
    
    D = spm_eeg_montage(S);
    
    % If you included bad channel interpolation, unset bad chans:
    if any(bcmontagefname)
        msg = 'and badchan-interpolated';
        fprintf('++ Unsetting interpolated bad channels\n');
        D = badchannels(D,D.badchannels,0);
        D.save;
    end

    fprintf('++ Cleaned %s data saved to %s\n',msg,D.fname);

end

%% Save in ICA:

ICA.removed = comp2rem;
ICA.montagefname = montagefname;
save(ICA.fname,'ICA');
fprintf('++ Output struct ICA saved to %s\n',ICA.fname);

return

