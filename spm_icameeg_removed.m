function [removed,report,fmain,fchild] = spm_icameeg_removed(S)
% Report removed ICA components and summarise evidence. 
%  FORMAT: [removed, report, fmain, fchild] = spm_icameeg_removed(S)
%  INPUT: Struct 'S' with fields:
%   S.D            - MEEG object or filename of MEEG object
%   S.ICA          - ICA struct (from spm_eeglab_runica)
%   S.dotopo       - plot topographies of removed components (def:1)
%   S.savetopo     - save topo figure as fig and png? (1(/prefix)|0) (def:'topo_removed')
%   S.doeegplot    - plot time-courses of removed components (def:0)
%   S.child.*      - see spm_eeglab_eegplot child plot options (def: EOG)
%
%  OUTPUT: 
%   removed        - removed artefact components (indices)
%   report         - table summarising removed components and evidence
%   fmain          - figure handle of main eegplot
%   fchild         - figure handle of child eegplot
%
%  spm_icameeg and spm_eeglab tools
%  by Jason Taylor (05/Aug/2021) jason.taylor@manchester.ac.uk

%--------------------------------------------------

% Check inputs:
if ~isstruct(S.ICA)
    load(S.ICA);
else
    ICA = S.ICA;
    S.ICA = ICA.fname; % save RAM
end
try dotopo = S.dotopo;       catch, dotopo = 1;    end
try savetopo = S.savetopo;   catch, savetopo = 'topo_removed';  end
try doeegplot = S.doeegplot; catch, doeegplot = 0; end
try childfield = S.child;    catch, childfield = []; end

% Load SPM-format data file:
D = spm_eeg_load(S.D);
[~,fstem] = fileparts(D.fname);

fprintf('\n\n');
fprintf('++ %s\n',datestr(now));
fprintf('++ RUNNING spm_icameeg_removed ON %s\n',ICA.fname);


%% Create report

removed = ICA.removed;

icamethods = fieldnames(ICA.artefact);
report = cell(0,5);

for r=1:length(removed)
    comp = removed(r);

    for m=1:length(icamethods)
        icamethod = icamethods{m};
    
        if strcmpi('manual',icamethod)
            ind = find([ICA.artefact.manual.suspect] == comp);
            if any(ind)
                reason = ICA.artefact.manual(ind).reason;
                report(end+1,:) = {comp 'manual' reason [] []};
            end
            
        else
        
            arts = fieldnames(ICA.artefact.(icamethod));
            
            for j=1:length(arts)
                art = arts{j};
                
                rmat = ICA.artefact.(icamethod).(art).summary;
                ind = find(rmat(:,1)==comp);
                r = rmat(ind,2);
                z = rmat(ind,3);
                report(end+1,:) = {comp icamethod art r z};
            end
 
        end
    end
end

report = sortrows(report,1);

% Save report:
reportfname = sprintf('report_removed_%s.mat',fstem);
save(reportfname,'report');
fprintf('++ Report saved to %s\n',reportfname);


%% Print report to screen

fprintf('ICA removed components report for %s\n',D.fname);

%removed = unique([report{:,1}]);
for comp=removed
    fprintf('IC%02d\n',comp);
    inds = find([report{:,1}] == comp);
    for i=inds
        fprintf('\tr=%.4f (z=%.4f) %s: %s\n',report{i,4},report{i,5},report{i,2},report{i,3});
    end
end

%% Plot component topographies and time-courses:

if dotopo
    
    fprintf('++ Plotting topographies of removed components\n')
    
    S=[];
    S.D              = D.fname;
    S.ICA            = ICA.fname;
    S.components     = removed;
    S.savefig        = savetopo;
    spm_icameeg_topo(S);
    
end

if doeegplot
    
    fprintf('++ Plotting time-courses of removed components and art chans\n');
    
    S=[];
    S.D              = D.fname;
    S.ICA            = ICA.fname;
    S.compinds       = removed;
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
    
