% Effects: Dilution, Group Detoxification, Group Sellout, Switchover, Best Case
% Worst Case

root_pathname = 'C:\Users\Wim\Documents\KE\Bsc Thesis\Code\Aggregation\Cluster Data\';
%folder_names = {'best_case', 'danger in numbers', 'dilution', 'group detoxification',...
%   'group sellout', 'switchover', 'worst_case'};
folder_names = {'danger in numbers'};
filenames = {'N=1','N=5','N=10'};

tmax = 100;

red = [1 0 0];
green = [0 1 0];
blue = [0 0 1];

for i=1:numel(folder_names)
    fol_name = sprintf('%s%s',root_pathname,folder_names{i});
    for j=1:numel(filenames)
        myfig=figure;
        title('Opt (Left) vs Cstnt (Right), Dilution Effect');
        trunk_filename = sprintf('%s%s\\alternate_setup\\%s_%s',root_pathname,folder_names{i},folder_names{i},filenames{j});
        opt_filename = sprintf('%s_DIFFERENT_ALPHA.m',trunk_filename);
        cons_filename = sprintf('%s_CONS_DIFFERENT_ALPHA.m',trunk_filename);
        opt=dlmread(opt_filename);
        cons=dlmread(cons_filename);

        muXopt=mean(opt(:,2)); % Mean population density
        muUopt=mean(opt(:,3)); % Mean resistance strategy

        x_label = sprintf('mu_{X} = %.3f',muXopt); % Label for population mean
        u_label = sprintf('mu_{u} = %.3f',muUopt); % Label for resistance mean

        % Optimized treatment results

        % Population subplot
        subplot(321),plot(opt(:,1),opt(:,2),'r'), line([0 100], [muXopt muXopt], 'Color', 'k')
        title('Population vs time (Optimized treatment)'),axis([0 tmax 0 100])
        text(100-100*0.2,muXopt-10,x_label)

        % Optimized treatment regime subplot
        subplot(323), treat1=area(opt(:,1),opt(:,4)),set(treat1,'FaceColor',green),
        title('Optimized treatment schedule'),axis([0 tmax 0 0.5*1.25])

        % Resistance strategy subplot
        res_ax_max = max(max(opt(:,3)),max(cons(:,3)));
        res_ax_max = res_ax_max * 1.1;
        subplot(325), plot(opt(:,1),opt(:,3),'b'), line([0 tmax], [muUopt muUopt], 'Color', 'k')
        title('Evolved resistance vs time (Optimized treatment)'),axis([0 tmax 0 res_ax_max]),
        text(100-100*0.2,muUopt+max(opt(:,3))/10,u_label)


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Constant treatment results %%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        muXcons=mean(cons(:,2)); % Mean population density
        muUcons=mean(cons(:,3)); % Mean resistance strategy

        x_label = sprintf('mu_{X} = %.3f',muXcons); % Label for population mean
        u_label = sprintf('mu_{u} = %.3f',muUcons); % Label for resistance mean
        %x_title = sprintf('Pop dty vs time, %s, alpha=%.1f, beta=%.1f, effect=%s',currN.name,params.alpha,params.beta,eff.name);

        % Population subplot
        subplot(322),plot(cons(:,1),cons(:,2),'r'), line([0 tmax], [muXcons muXcons], 'Color', 'k')
        title('Population vs time (Constant treatment)'),axis([0 tmax 0 100])
        text(100-100*0.2,muXcons-10,x_label)

        % Optimized treatment regime subplot
        subplot(324), treat2=area(cons(:,1),cons(:,4)),set(treat2,'FaceColor',green),
        title('Constant treatment schedule'),axis([0 tmax 0 0.5*1.25])

        % Resistance strategy subplot
        subplot(326), plot(cons(:,1),cons(:,3),'b'), line([0 tmax], [muUcons muUcons], 'Color', 'k')
        title('Evolved resistance vs time (Constant treatment)'),axis([0 tmax 0 res_ax_max]),
        text(100-100*0.2,muUcons+max(cons(:,3))/10,u_label)

        
        new_filename=strrep(sprintf('%s %s',folder_names{i},filenames{j}),'_',' ');
        p=mtit(new_filename);
    end
end