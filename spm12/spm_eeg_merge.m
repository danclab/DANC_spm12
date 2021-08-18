function Dout = spm_eeg_merge(S)
% Concatenate epoched single trial files
% FORMAT Dout = spm_eeg_merge(S)
%
% S           - input structure (optional)
%  fields of S:
%   S.D       - character array containing filename of M/EEG mat-files
%               or cell array of D's
%   S.recode  - this field specifies how the condition labels will be
%               translated from the original files to the merged file.
%               Several options are possible:
%                 'same'        - leave the condition labels unchanged
%                 'addfilename' - add the original file name to condition
%                                 label
%                 old way specification - (for backward compatibility)                     
%                       a cell array where each cell contains a condition
%                       label. The ordering of these labels must be such 
%                       that each row in the cell matrix specifies the 
%                       conditionlabels for one of the selected files.
%                 specification via recoding rules - for this S.recode
%                       should be a structure array where each element 
%                       specifies a rule using the following fields:
%                            file - can be a cell array of strings with 
%                                   file names, a vector of file indices 
%                                   or a string with regular expression 
%                                   matching the files to which the rule 
%                                   will apply.
%                            labelorg - can be a cell array of condition 
%                                   labels or a string with regular 
%                                   expression matching the condition 
%                                   labels to which this rule will apply.
%                            labelnew - new label for the merged file. It
%                                   can contain special tokens #file# and
%                                   #labelorg# that will be replaced by 
%                                   the original file name and original 
%                                   condition label respectively.
%                       The rule will be applied one after the other so 
%                       the last rule takes precedences. Trials not 
%                       matched by any of the rules will keep their 
%                       original labels.
%                       Example:
%                          S.recode(1).file     = '.*';
%                          S.recode(1).labelorg = '.*';
%                          S.recode(1).labelnew = '#labelorg# #file#';
%                       has the same effect as the 'addfilename' option.
%   S.prefix  - prefix for the output file (default - 'c')
%
% 
% Dout        - MEEG object (also written to disk)
%__________________________________________________________________________
%
% This function can be used to merge M/EEG files to one file. This is
% useful whenever the data are distributed over multiple files, but one
% wants to use all information in one file. For example, when displaying
% data (SPM displays data from only one file at a time), or merging
% information that has been measured in multiple sessions.
%__________________________________________________________________________
% Copyright (C) 2008-2017 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel, Vladimir Litvak, Doris Eckstein, Rik Henson
% $Id: spm_eeg_merge.m 7125 2017-06-23 09:49:29Z guillaume $

SVNrev = '$Rev: 7125 $';

%-Startup
%--------------------------------------------------------------------------
spm('FnBanner', mfilename, SVNrev);
spm('FigName','M/EEG Merge'); spm('Pointer','Watch');

if ~isfield(S,'prefix'), S.prefix = 'c';    end
if ~isfield(S,'recode'), S.recode = 'same'; end

%-Load MEEG data
%--------------------------------------------------------------------------
D = S.D;

if ischar(D)
    F = cell(1,size(D,1));
    try
        for i = 1:size(D, 1)
            F{i} = spm_eeg_load(deblank(D(i, :)));
        end
        D = F;
    catch exception
        disp(getReport(exception));
        error('Trouble reading files.');
    end
end

Nfiles = length(D);

if Nfiles < 2
    %error('Need at least two files for merging.');
end

%-Check input and determine number of new number of trial types
%--------------------------------------------------------------------------
Ntrials = [];
megsens = [];
eegsens = [];
fid     = [];
isTF    =  strncmpi(D{1}.transformtype,'TF',2); % TF and TFphase

for i = 1:Nfiles
    if ~isequal(D{i}.transformtype, D{1}.transformtype)
        error(['The datasets do not contain the same kind of data.\n'...
               'There is a difference between files\n\t%s\nand\n\t%s.'], ...
               D{1}.fname, D{i}.fname);
    end

%     if D{1}.nchannels ~= D{i}.nchannels
%         error(['Data don''t have the same number of channels.\n' ...
%                'There is a difference between files\n\t%s\nand\n\t%s.'], ...
%                D{1}.fname, D{i}.fname);
%     end

    if D{1}.nsamples ~= D{i}.nsamples
        error(['Data don''t have the same number of time points.\n' ...
               'There is a difference between files\n\t%s\nand\n\t%s.'], ...
               D{1}.fname, D{i}.fname);
    end

    if D{1}.fsample ~= D{i}.fsample
        error(['Data don''t have the same sampling rate.\n' ...
               'There is a difference between files\n\t%s\nand\n\t%s.'], ...
               D{1}.fname, D{i}.fname);
    end

    if isTF &&  ~isequal(D{1}.frequencies, D{i}.frequencies)
        error(['Data don''t have the same frequencies.\n' ...
               'There is a difference between files\n\t%s\nand\n\t%s.'], ...
               D{1}.fname, D{i}.fname);
    end

    if ~isempty(D{i}.sensors('MEG'))
        megsens = spm_cat_struct(megsens, D{i}.sensors('MEG'));
    end
    
    if ~isempty(D{i}.sensors('EEG'))
        eegsens = spm_cat_struct(eegsens, D{i}.sensors('EEG'));
    end
    
    if ~isempty(megsens) || ~isempty(eegsens)
        fid = spm_cat_struct(fid, D{i}.fiducials);
    end
    
    Ntrials = [Ntrials D{i}.ntrials];
end

%-Prepare some useful lists
%--------------------------------------------------------------------------
F        = {};
Find     = [];
clb      = {};
for i = 1:Nfiles
    F{i} = fname(D{i});
    clb  = [clb D{i}.conditions];
    Find = [Find i*ones(1, D{i}.ntrials)];
end
uclb     = unique(clb);


%-Specify condition labels recoding
%--------------------------------------------------------------------------
if ~isfield(S, 'recode')
    S.recode = spm_input('What to do with condition labels?', 1, 'm',...
        'Leave as they are|Add file name|Specify rules for recoding|Specify recoding the old way', strvcat('same', 'addfilename', 'rules', 'old'));
end

if isequal(S.recode, 'old')
    S.recode = {};
    for i = 1:Nfiles
        for j = 1:nconditions(D{i})
            S.recode{i}{j} = spm_input(sprintf('Labels: %s', spm_file(D{i}.fname, 'basename')),...
                '+1', 's', D{i}.condlist{j});
        end
    end
elseif isequal(S.recode, 'rules')
    S.recode = [];
    stop     = 0;
    ind      = 1;   
    while ~stop
        spm_input(['Please define rule ' num2str(ind) ':'], 1, 'd');
        switch spm_input('To which files will this rule apply?', '+1', 'm',...
                'All the files|Specify indices|Select files|Wildcard expression (*,?)|Regular expression',...
                strvcat('all', 'indices', 'select', 'wildcard', 'regexp'))
            case 'all'
                S.recode(ind).file = '.*';
            case 'indices'
                S.recode(ind).file =  spm_input('Input file indices', '+1', 'n', num2str(ind), [1 Inf]);
            case 'select'
                [selection, ok]= listdlg('ListString', F, 'SelectionMode', 'multiple' ,'Name', 'Select files' , 'ListSize', [400 300]);
                if ok
                    S.recode(ind).file = F(selection);
                else
                    continue;
                end
            case 'wildcard'
                S.recode(ind).file = regexptranslate('wildcard' , spm_input('Input wildcard expresssion', '+1', 's',  '*'));
            case 'regexp'
                S.recode(ind).file = spm_input('Input regular expresssion', '+1', 's',  '.*');
        end

        switch spm_input('What conditions will be renamed?', '+1', 'm',...
                'All|Select|Specify by wildcard expression (*,?)|Specify by regular expression',...
                strvcat('all', 'select', 'wildcard', 'regexp'))
            case 'all'
                S.recode(ind).labelorg = '.*';
            case 'select'
                [selection, ok]= listdlg('ListString', uclb, 'SelectionMode', 'multiple' ,'Name', 'Select conditions' , 'ListSize', [400 300]);
                if ok
                    S.recode(ind).labelorg = uclb(selection);
                else
                    continue;
                end
            case 'wildcard'
                S.recode(ind).labelorg = regexptranslate('wildcard' , spm_input('Input wildcard expresssion', '+1', 's',  '*'));
            case 'regexp'
                S.recode(ind).labelorg = spm_input('Input regular expresssion', '+1', 's',  '.*');
        end

        S.recode(ind).labelnew = spm_input('Input the new name?', '+1', 's',  '');
     
        stop = spm_input('Define another rule?','+1','yes|stop', [0 1], 0);
        ind  = ind+1;
    end
end

%-Generate new meeg object with new filename
%--------------------------------------------------------------------------
Dout = D{end};
newfilename = spm_file(fnamedat(Dout), 'path',pwd, 'prefix',S.prefix);

if ~isTF
    Dout = clone(Dout, newfilename, [D{end}.nchannels Dout.nsamples sum(Ntrials)]);
else
    Dout = clone(Dout, newfilename, [D{end}.nchannels Dout.nfrequencies Dout.nsamples sum(Ntrials)]);
end


%-Perform condition labels recoding
%--------------------------------------------------------------------------
if isequal(S.recode, 'same')
    Dout = conditions(Dout, ':', clb);
elseif isequal(S.recode, 'addfilename')
    for i = 1:numel(clb)
        clb{i} = [clb{i} ' ' spm_file(F{Find(i)}, 'basename')];
    end
    Dout = conditions(Dout, ':', clb);
elseif iscell(S.recode)
    for i = 1:Nfiles
        ind = find(Find == i);        
               
        for j = 1:D{i}.nconditions
            clb(ind(strmatch(D{i}.condlist{j}, clb(ind), 'exact'))) = S.recode{i}(j);
        end
    end
    Dout = conditions(Dout, ':', clb);
elseif isstruct(S.recode)
    clbnew = clb;
    
    for i = 1:numel(S.recode)
        if isnumeric(S.recode(i).file)
            ind = S.recode(i).file;
        elseif iscell(S.recode(i).file)
            ind = spm_match_str(F, S.recode(i).file);
        elseif ischar(S.recode(i).file)
            ind = find(~cellfun('isempty', regexp(F, S.recode(i).file)));
        else
            error('Invalid file specification in recoding rule.');
        end
        
        ind = find(ismember(Find, ind));
        
        if iscell(S.recode(i).labelorg)
            ind = ind(ismember(clb(ind), S.recode(i).labelorg));
        elseif ischar(S.recode(i).labelorg)
            ind = ind(~cellfun('isempty', regexp(clb(ind), S.recode(i).labelorg)));
        else
            error('Invalid original condition label specification in recoding rule.');
        end
        
        for j = 1:length(ind)
            labelnew = S.recode(i).labelnew;
            labelnew = strrep(labelnew, '#file#', spm_file(F{Find(ind(j))}, 'basename'));
            labelnew = strrep(labelnew, '#labelorg#', clb(ind(j)));
            
            clbnew{ind(j)} = labelnew;
        end
    end
    
    Dout = conditions(Dout, ':', clbnew);
end
            
%-Average sensor locations
%--------------------------------------------------------------------------
CmdLine = spm('CmdLine');
if CmdLine, h = []; end
if ~isempty(megsens)
    if ~CmdLine, spm_figure('GetWin','Graphics');clf; end
    if ~isempty(eegsens)
        if ~CmdLine, h = subplot(2, 1, 1); end
        aeegsens = ft_average_sens(eegsens, 'weights', Ntrials, 'feedback', h);
        Dout = sensors(Dout, 'EEG', aeegsens);
        
        if ~CmdLine, h = subplot(2, 1, 2); end
    else
        if ~CmdLine, h = axes; end
    end
    
    common_chans=intersect(megsens(1).label, megsens(2).label);
    for i=3:length(megsens)
        common_chans=intersect(common_chans,megsens(i).label);
    end
    min_coils=Inf;
    min_coil_idx=-1;
    for i=1:length(megsens)
        if length(megsens(i).label)<min_coils
            min_coils=length(megsens(i).label);
            min_coil_idx=i;
        end
    end
    

    for idx=1:length(megsens)
        chans_to_remove=setdiff(megsens(idx).label,common_chans);

        for j=1:length(chans_to_remove)
            chan_to_remove=chans_to_remove{j};

            nchans=size(megsens(idx).chanori,1);
            chan=find(strcmp(megsens(idx).label,chan_to_remove));
            chan_idx=setdiff([1:nchans],[chan]);
            megsens(idx).chanori=megsens(idx).chanori(chan_idx,:);
            megsens(idx).chanpos=megsens(idx).chanpos(chan_idx,:);
            megsens(idx).chantype=megsens(idx).chantype(chan_idx);
            megsens(idx).chantypeold=megsens(idx).chantypeold(chan_idx);
            megsens(idx).chanunit=megsens(idx).chanunit(chan_idx);
            megsens(idx).chanunitold=megsens(idx).chanunitold(chan_idx);
            megsens(idx).label=megsens(idx).label(chan_idx);
            megsens(idx).labelold=megsens(idx).labelold(chan_idx);
            %megsens(idx).tra=megsens(idx).tra(chan_idx,:);

            nchans=length(megsens(idx).balance.G1BR.labelold);
            chan=find(strcmp(megsens(idx).balance.G1BR.labelold,chan_to_remove));
            chan_idx=setdiff([1:nchans],[chan]);
            megsens(idx).balance.G1BR.labelold=megsens(idx).balance.G1BR.labelold(chan_idx);
            megsens(idx).balance.G1BR.labelnew=megsens(idx).balance.G1BR.labelnew(chan_idx);
            megsens(idx).balance.G1BR.tra=megsens(idx).balance.G1BR.tra(chan_idx,chan_idx);
            megsens(idx).balance.G1BR.chantypeold=megsens(idx).balance.G1BR.chantypeold(chan_idx);
            megsens(idx).balance.G1BR.chanunitold=megsens(idx).balance.G1BR.chanunitold(chan_idx);
            megsens(idx).balance.G1BR.chantypenew=megsens(idx).balance.G1BR.chantypenew(chan_idx);
            megsens(idx).balance.G1BR.chanunitnew=megsens(idx).balance.G1BR.chanunitnew(chan_idx);

            nchans=length(megsens(idx).balance.G2BR.labelold);
            chan=find(strcmp(megsens(idx).balance.G2BR.labelnew,chan_to_remove));
            chan_idx=setdiff([1:nchans],[chan]);
            megsens(idx).balance.G2BR.labelold=megsens(idx).balance.G2BR.labelold(chan_idx);
            megsens(idx).balance.G2BR.labelnew=megsens(idx).balance.G2BR.labelnew(chan_idx);
            megsens(idx).balance.G2BR.tra=megsens(idx).balance.G2BR.tra(chan_idx,chan_idx);
            megsens(idx).balance.G2BR.chantypeold=megsens(idx).balance.G2BR.chantypeold(chan_idx);
            megsens(idx).balance.G2BR.chanunitold=megsens(idx).balance.G2BR.chanunitold(chan_idx);
            megsens(idx).balance.G2BR.chantypenew=megsens(idx).balance.G2BR.chantypenew(chan_idx);
            megsens(idx).balance.G2BR.chanunitnew=megsens(idx).balance.G2BR.chanunitnew(chan_idx);

            nchans=length(megsens(idx).balance.G3BR.labelold);
            chan=find(strcmp(megsens(idx).balance.G3BR.labelold,chan_to_remove));
            chan_idx=setdiff([1:nchans],[chan]);
            megsens(idx).balance.G3BR.labelold=megsens(idx).balance.G3BR.labelold(chan_idx);
            megsens(idx).balance.G3BR.labelnew=megsens(idx).balance.G3BR.labelnew(chan_idx);
            megsens(idx).balance.G3BR.tra=megsens(idx).balance.G3BR.tra(chan_idx,chan_idx);
            megsens(idx).balance.G3BR.chantypeold=megsens(idx).balance.G3BR.chantypeold(chan_idx);
            megsens(idx).balance.G3BR.chanunitold=megsens(idx).balance.G3BR.chanunitold(chan_idx);
            megsens(idx).balance.G3BR.chantypenew=megsens(idx).balance.G3BR.chantypenew(chan_idx);
            megsens(idx).balance.G3BR.chanunitnew=megsens(idx).balance.G3BR.chanunitnew(chan_idx);

            noldchans=length(megsens(idx).balance.custom.labelold);
            nnewchans=length(megsens(idx).balance.custom.labelnew);
            oldchan=find(strcmp(megsens(idx).balance.custom.labelold,chan_to_remove));
            newchan=find(strcmp(megsens(idx).balance.custom.labelnew,chan_to_remove));
            old_chan_idx=setdiff([1:noldchans],[oldchan]);
            new_chan_idx=setdiff([1:nnewchans],[newchan]);
            megsens(idx).balance.custom.labelold=megsens(idx).balance.custom.labelold(old_chan_idx);
            megsens(idx).balance.custom.labelnew=megsens(idx).balance.custom.labelnew(new_chan_idx);
            megsens(idx).balance.custom.tra=megsens(idx).balance.custom.tra(new_chan_idx,old_chan_idx);
            megsens(idx).balance.custom.chantypeold=megsens(idx).balance.custom.chantypeold(old_chan_idx);
            megsens(idx).balance.custom.chanunitold=megsens(idx).balance.custom.chanunitold(old_chan_idx);
            megsens(idx).balance.custom.chantypenew=megsens(idx).balance.custom.chantypenew(new_chan_idx);
            megsens(idx).balance.custom.chanunitnew=megsens(idx).balance.custom.chanunitnew(new_chan_idx);

        end
        megsens(idx).coilpos=megsens(min_coil_idx).coilpos;
        megsens(idx).coilori=megsens(min_coil_idx).coilori;
        megsens(idx).tra=megsens(min_coil_idx).tra;
    end

    [amegsens,afid] = ft_average_sens(megsens, 'fiducials', fid, 'weights', Ntrials, 'feedback', h);
    Dout = sensors(Dout, 'MEG', amegsens);
    Dout = fiducials(Dout, afid);
elseif ~isempty(eegsens)
    if ~CmdLine, spm_figure('GetWin','Graphics');clf; end
    if ~CmdLine, h = axes; end
    [aeegsens,afid] = ft_average_sens(eegsens, 'fiducials', fid, 'weights', Ntrials, 'feedback', h);
    Dout = sensors(Dout, 'EEG', aeegsens);
    Dout = fiducials(Dout, afid);
end

%-Write files
%--------------------------------------------------------------------------
spm_progress_bar('Init', Nfiles, 'Files merged');

k = 0;

for i = 1:Nfiles

    ind = union(Dout.badchannels, D{i}.badchannels);
    if ~isempty(ind)
        Dout = badchannels(Dout, ind, 1);
    end

    chan_inds=[];
    for j=1:length(Dout.chanlabels)
        chan_idx=find(strcmp(D{i}.chanlabels,Dout.chanlabels(j)));
        if length(chan_idx)
            chan_inds(end+1)=chan_idx;
        end
    end
    % write trial-wise to avoid memory mapping error
    for j = 1:D{i}.ntrials
        k = k + 1;
        if ~isTF
            %Dout(1:Dout.nchannels, 1:Dout.nsamples, k) =  D{i}(chan_inds,:,j);
            for l=1:length(Dout.chanlabels)
               chan_idx=find(strcmp(D{i}.chanlabels,Dout.chanlabels(l)));
               if length(chan_idx)
                   Dout(l, 1:Dout.nsamples, k) =  D{i}(chan_idx,:,j);
               end
            end
        else
            Dout(1:Dout.nchannels, 1:Dout.nfrequencies, 1:Dout.nsamples, k) =  D{i}(chan_inds,:,:,j);
        end
        Dout = badtrials(Dout, k, badtrials(D{i}, j));
    end
    
    % Propagate some useful information from the original files to the
    % merged file
    Dout = repl(Dout, find(Find == i), D{i}.repl);
    Dout = trialonset(Dout, find(Find == i), D{i}.trialonset);
    Dout = trialtag(Dout, find(Find == i), D{i}.trialtag);
    Dout = events(Dout, find(Find == i), D{i}.events);
    
    spm_progress_bar('Set', i);

end

%-Save new M/EEG data
%--------------------------------------------------------------------------
Dout = Dout.history('spm_eeg_merge', S, 'reset');
save(Dout);

%-Cleanup
%--------------------------------------------------------------------------
spm_progress_bar('Clear');
fprintf('%-40s: %30s\n','Completed',spm('time'));                       %-#
spm('FigName','M/EEG merge: done'); spm('Pointer','Arrow');
