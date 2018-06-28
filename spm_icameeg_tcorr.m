function ICA = spm_icameeg_tcorr(S)
% Correlate ICA time-courses with artefact channels.
%  FORMAT: ICA = spm_icameeg_tcorr(S)
%  INPUT: Struct 'S' with fields:
%   S.D          - MEEG object or filename of MEEG object
%   S.ICA        - ICA struct or filename (output of spm_eeglab_runica)
%   S.artchans   - Artefact chan names, cell array (default: {'VEOG'})
%   S.artfilts   - Bandpass filters to use (default: {[1 20]};
%   S.thresh     - Threshold for z-score of temporal correlation (def:2)
%  OUTPUT: 
%   ICA          - ICA struct with 'artefact.tcorr' field with subfields:
%    (artchan).summary - correlation matrix with columns:
%                         component index, r, z-score of r
%    (artchan).suspect - indices of suprathreshold components
%    (artchan).thresh  - threshold used
%
%  spm_icameeg and spm_eeglab tools
%  by Jason Taylor (09/Mar/2018) jason.taylor@manchester.ac.uk

%  NOTE: need to sort out topo of MEGPLANAR sensors (RMS at each location)

%--------------------------------------------------

%% Check inputs:

if ~isstruct(S.ICA)
    load(S.ICA);
else
    ICA = S.ICA;
    S.ICA = ICA.fname; % save RAM
end
try artchans  = S.artchans;  catch, artchans  = {'VEOG'}; end
try artfilts  = S.artfilts;  catch, artfilts  = {[1 20]}; end
try thresh    = S.thresh;    catch, thresh    = 2;        end

%% Load SPM-format data file:

D = spm_eeg_load(S.D);
[~,fstem] = fileparts(D.fname);

fprintf('\n\n');
fprintf('++ %s\n',datestr(now));
fprintf('++ RUNNING spm_icameeg_tcorr ON %s\n',D.fname);
fprintf('++ USING: thresh =  %g\n',thresh);
artchantxt = sprintf('%s ',artchans{:});
fprintf('++ USING: artchans: %s\n',artchantxt);
artfilttxt = sprintf('%g-%g ',artfilts{:});
fprintf('++ USING: filters: %s\n',artfilttxt);

% Prepare figure:
spm_figure('Clear','Graphics');
fig = spm_figure('GetWin','Graphics');

for i=1:length(artchans)
    artchan = artchans{i};
    artfilt = artfilts{i};
    suspect = [];
    
    %% Extract data, filter it:
    if ~iscell(ICA.timewin)
        a = selectdata(D,artchan,ICA.timewin,[]);
    else
        a = [];
        for tw=1:length(ICA.timewin)
            a = [a selectdata(D,artchan,ICA.timewin{tw},[])];
        end
    end
    if ~strcmpi(D.type,'continuous')
        % NOTE: It may be dodgy to filter concatenated epochs!
        a = reshape(a,size(a,1),size(a,2)*size(a,3));
    end
    af = ft_preproc_bandpassfilter(a,D.fsample,artfilt,[],'but','twopass','reduce');
    
    % Filter component activations:
    actf = ft_preproc_bandpassfilter(ICA.activations(:,:),D.fsample,artfilt,[],'but','twopass','reduce'); 

    %% Compute temporal correlations:
    rt = zeros(ICA.ncomp,1);
    for j=1:ICA.ncomp
        try
            rt(j) = abs(corr(af',actf(j,:)'));
        catch
            tmp = corrcoef(af',actf(j,:)');
            rt(j) = abs(tmp(2));
        end
    end
    
    % Compute z-score of correlations (for thresholding)
    zrt = (rt-mean(rt))/std(rt);

    % Construct summary matrix, sort by z-score (descending):
    rmat = [(1:size(actf,1))' rt zrt];
    rmat = sortrows(rmat,2);
    
    %% Plot z-scores:
    subplot(length(artchans),1,i); 
    bar(zrt); title(artchan);
    xlabel('IC'); ylabel('z-score of correlation');
    set(gca,'xtick',1:ICA.ncomp)
    set(gca,'fontsize',8)
    axis tight

    %% Report suprathreshold correlations:
    indsupthresh = find(rmat(:,3)>thresh);
    fprintf('++ Found %d suprathreshold temporal correlations with %s\n',length(indsupthresh),artchan);
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

    %% Store it:
    ICA.artefact.tcorr.(artchan).summary = rmat;
    ICA.artefact.tcorr.(artchan).suspect = suspect;
    ICA.artefact.tcorr.(artchan).thresh  = thresh;
    
end

%% Save ICA:
save(ICA.fname,'ICA');
fprintf('++ Saved output ICA struct to %s\n',ICA.fname);

% Save z-score figure:
artchantxt = sprintf('%s_',artchans{:});
figfname = sprintf('summary_tcorr_%s%s.png',artchantxt,fstem);
print(fig,'-dpng',figfname);
fprintf('++ Saved figure to %s\n',figfname);

return
