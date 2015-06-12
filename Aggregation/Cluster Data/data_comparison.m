% for every folder i
%   for every size of N j in folder i
%       read in optimized file o
%       read in constant file c
%       compare o with c
%       write to cell array of structs, one struct for every effect

%folder_names = {'best_case', 'danger in numbers', 'dilution', 'group detoxification',...
%    'group sellout', 'switchover', 'worst_case'};
clear variables

root_name = 'C:\Users\Wim\Documents\KE\Bsc Thesis\Code\Aggregation\Cluster Data\';
folder_names = {'danger in numbers'};
filenames = {'N=1', 'N=5', 'N=10'};
rep_filenames = {};
for i=1:numel(filenames)
    rep_filenames{i} = strrep(filenames{i},'=','');
end

results = {};
for i=1:numel(folder_names)
    curr_folder_name = folder_names{i};
    %res = struct(rep_filenames{1},[],rep_filenames{2},[],rep_filenames{3},[]);
    effect_res = struct('name',curr_folder_name);
    trunk_filename = sprintf('%s%s\\alternate_setup\\%s',root_name,curr_folder_name,curr_folder_name);
    for j=1:numel(filenames)
        curr_file_name = filenames{j};
        opt_filename = sprintf('%s_%s.m',trunk_filename,curr_file_name);
        opt = dlmread(opt_filename);
        
        cons_filename = sprintf('%s_%s_CONS.m',trunk_filename,curr_file_name);
        cons = dlmread(cons_filename);
        
        muXopt = mean(opt(:,2));
        muXcons = mean(cons(:,2));
        
        muUopt = mean(opt(:,3));
        muUcons = mean(cons(:,3));
        
        N_res = struct('name', rep_filenames{j}, 'muXopt', muXopt, 'muXcons', muXcons, 'diff', 100*((muXcons-muXopt)/muXcons));
        effect_res = setfield(effect_res, rep_filenames{j}, N_res);
    end
    results{i} = effect_res;
end