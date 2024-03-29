function [suspects,report,fmain,fchild] = spm_icameeg_suspects(S)
% Report 'suspect' ICA components from various methods. 
%  FORMAT: [suspects, report] = spm_icameeg_suspects(S)
%  INPUT: Struct 'S' with fields:
%   S.D            - MEEG object or filename of MEEG object
%   S.ICA          - ICA struct (from spm_eeglab_runica)
%   S.dotopo       - plot topographies of suspect components (def:0)
%   S.doeegplot    - plot time-courses of suspect components (def:0)
%   S.child.*      - see spm_eeglab_eegplot child plot options (def: EOG)
%
%  OUTPUT: 
%   suspects       - suspected artefact components
%   report         - table summarising suspects and reasons
%   fmain          - figure handle of main eegplot
%   fchild         - figure handle of child eegplot
%
%  spm_icameeg and spm_eeglab tools
%  by Jason Taylor (09/Mar/2018) jason.taylor@manchester.ac.uk
%   + (07/Jun/2021) jt: added eegplot child options, fmain fchild output

%--------------------------------------------------

% Check inputs:
if ~isstruct(S.ICA)
    load(S.ICA);
else
    ICA = S.ICA;
    S.ICA = ICA.fname; % save RAM
end
try dotopo = S.dotopo;       catch, dotopo = 0;    end
try doeegplot = S.doeegplot; catch, doeegplot = 0; end
try childfield = S.child;    catch, childfield = []; end

% Load SPM-format data file:
D = spm_eeg_load(S.D);
[~,fstem] = fileparts(D.fname);

fprintf('\n\n');
fprintf('++ %s\n',datestr(now));
fprintf('++ RUNNING spm_icameeg_suspects ON %s\n',ICA.fname);


%% Create report

icamethods = fieldnames(ICA.artefact);
report = cell(0,5);
for m=1:length(icamethods)
    icamethod = icamethods{m};
    
    if strcmpi('manual',icamethod)
        for j=1:length(ICA.artefact.manual)
            if ~isempty(ICA.artefact.manual(j).suspect)
                comp = ICA.artefact.manual(j).suspect;
                reason = ICA.artefact.manual(j).reason;
                report(end+1,:) = {comp 'manual' reason [] []};
            end
        end
    else
        
        arts = fieldnames(ICA.artefact.(icamethod));
        for j=1:length(arts)
            art = arts{j};
            suspect = ICA.artefact.(icamethod).(art).suspect;
            if any(suspect)
                rmat = sortrows(ICA.artefact.(icamethod).(art).summary,1);
                for j=1:length(suspect)
                    comp = suspect(j);
                    r = rmat(comp,2);
                    z = rmat(comp,3);
                    report(end+1,:) = {comp icamethod art r z};
                end
            end
        end
    end
end
report = sortrows(report,1);

% Save report:
reportfname = sprintf('report_suspects_%s.mat',fstem);
save(reportfname,'report');
fprintf('++ Report saved to %s\n',reportfname);


%% Print report to screen

fprintf('ICA ''suspect'' components report for %s\n',D.fname);

suspects = unique([report{:,1}]);
for comp=suspects
    fprintf('IC%02d\n',comp);
    inds = find([report{:,1}] == comp);
    for i=inds
        fprintf('\tr=%.4f (z=%.4f) %s: %s\n',report{i,4},report{i,5},report{i,2},report{i,3});
    end
end

%% Plot component topographies and time-courses:

if dotopo
    
    fprintf('++ Plotting topographies of suspect components\n')
    
    S=[];
    S.D              = D.fname;
    S.ICA            = ICA.fname;
    S.components     = suspects;
    S.savefig        = 0;
    spm_icameeg_topo(S);
    
end

if doeegplot
    
    fprintf('++ Plotting time-courses of suspect components and art chans\n');
    
    S=[];
    S.D              = D.fname;
    S.ICA            = ICA.fname;
    S.compinds       = suspects;
    S.title          = 'IC Activations';
    S.spacing        = 10;
    if ~isempty(childfield)
        S.child = childfield;
    else
        S.child.chantype = {'EOG'};
        S.child.title    = 'EOG';
        S.child.spacing  = 300;
    end
    [fmain,fchild] = spm_eeglab_eegplot(S);
    
end

return
    
