spm_icameeg
spm_eeglab
spm_uomeeg
Toolboxes (if we can call them that?)
(09/Mar/2018) by Jason Taylor jason.taylor@manchester.ac.uk
(02/Jun/2018) JT: some minor bugfixes, check MEG compatibility
(25/Jun/2018) JT: name changes
 
For example usage, see:
script_spm_icameeg_long.m (best run cell-by-cell)
script_spm_icameeg_short.m (runnable start-to-finish as a script)

A very basic algorithm for removal of ICA components from SPM12 data using 
the montage function is here:
notes_basic_ica_montage_algorithm.m

To use spm_uomeeg_channelrepair.m you must copy ft_channelrepair_jt.m to 
<spm>/external/fieldtrip/
(where <spm> is the base spm directory, e.g., given by 'which spm' in matlab)

NOTE: spm_uomeeg_findbadchans.m is highly experimental at this point. Use with 
extreme caution! Also, it needs the stats toolbox :-\

NOTE: I think these functions now work with MEG data as well as EEG. 
However, topo plotting for MEGPLANAR (using RMS) has not been implemented.
