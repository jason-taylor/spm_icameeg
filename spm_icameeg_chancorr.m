function ICA = spm_icameeg_chancorr(S)
%  FORMAT: ICA = spm_icameeg_chancorr(S)
%  INPUT: Struct 'S' with fields:
%   S.D          - MEEG object or filename of MEEG object
%   S.ICA        - ICA struct or filename (output of spm_eeglab_runica)
%   S.thresh     - threshold for ratio max:next-highest corr (default: 10)
%  OUTPUT: 
%   ICA          - ICA struct: artefact.chancorr field with subfields:
%    (chantype).summary  - correlation matrix with columns (r) with
%                           channels in rows and components in columns.
%    (chantype).suspect  - indices of suprathreshold components
%    (chantype).thresh   - threshold used
%    (chantype).rmat     - correlation matrix
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
compinds = 1:ICA.ncomp;
try thresh = S.thresh; catch, thresh = 10; end

%% Load SPM-format data file, get channel info:
D = spm_eeg_load(S.D);
[~,fstem] = fileparts(D.fname);
nchan = length(ICA.chans);

fprintf('\n\n');
fprintf('++ %s\n',datestr(now));
fprintf('++ RUNNING spm_icameeg_chancorr ON %s\n',D.fname);
chantxt = sprintf('%s ',ICA.chans{:});
fprintf('++ USING: thresh =  %g\n',thresh);
fprintf('++ USING: channels: %s\n',chantxt);


%% Compute correlations between IC and channel time-courses:
fprintf('Progress (/%d)\n',ICA.ncomp);
rmat = zeros(nchan,ICA.ncomp);
for i=compinds
    fprintf('%d.',i);
    act = ICA.activations(i,:);
    for j=1:nchan

        % Extract data (check for timewin)
        if ~iscell(ICA.timewin)
            d = selectdata(D,ICA.chans{j},ICA.timewin,[]);
        else
            d = [];
            for tw=1:length(ICA.timewin)
                d = [d selectdata(D,ICA.chans{j},ICA.timewin{tw},[])];
            end
        end
        
        % Reshape to continuous if epoched:
        if D.ntrials>1
            d = reshape(d,1,size(d,2)*size(d,3));
        end
        
        % Correlations:
        try
            rmat(j,i) = abs(corr(d',act'));
        catch
            tmp = corrcoef(d',act');
            rmat(j,i) = abs(tmp(2));
        end
        
    end
end
fprintf('\n');

%% Plot correlations as scaled image:
%figure('color','w');
spm_figure('Clear','Graphics');
fig = spm_figure('GetWin','Graphics');
imagesc(rmat); set(gca,'clim',[0 1]); colormap jet
xlabel('IC'); ylabel('channel');
set(gca,'xtick',compinds); set(gca,'xticklabel',compinds);
set(gca,'ytick',(1:length(ICA.chans))); set(gca,'yticklabel',ICA.chans);
set(gca,'fontsize',8)
% Save figure as image:
figfname = sprintf('summary_chancorr_eeg_%s.png',fstem);
print(fig,'-dpng',figfname);
fprintf('++ Figure saved to %s\n',figfname);

%% Find components with a high ratio of max:next-highest chan corr
% (indicative of a component capturing noise on a single channel):
maxratio = [];
for i=compinds
    tmp = sort(rmat(:,i));
    maxratio(i) = tmp(end)/tmp(end-1);
end

% Report:
suspect = [];
suprathreshratio = find(maxratio>thresh);
if any(suprathreshratio)
    suspect(1,:) = suprathreshratio;
    fprintf('++ Found %d suprathreshold max:next-highest ratio of channel correlations\n',length(suprathreshratio));
    for i=suprathreshratio
        stchan = char(ICA.chans(find(rmat(:,i)==max(rmat(:,i)))));
        fprintf('IC%d \tratio=%0.1f\tr=%0.4f\tchan=%s\n',i,maxratio(i),max(rmat(:,i)),stchan);
    end
else
    fprintf('++ No suprathreshold max:next-highest ratio chan corr found.\n');
end

z = (maxratio-mean(maxratio))/std(maxratio)';
s = [(1:ICA.ncomp)' maxratio' z'];
s = sortrows(s,3);

%% Save in ICA:
ICA.artefact.chancorr.eeg.summary  = s;
ICA.artefact.chancorr.eeg.suspect  = suspect;
ICA.artefact.chancorr.eeg.thresh   = thresh;
ICA.artefact.chancorr.eeg.rmat     = rmat;

% Save ICA:
save(ICA.fname,'ICA');
fprintf('++ Saved output ICA struct to %s\n',ICA.fname);

return
