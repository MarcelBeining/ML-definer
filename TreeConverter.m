function [all_trees,MainDir] = TreeConverter(options)
if nargin < 1
    options = [];
end
MainDir = uigetdir([],'Select Root Directory of Animal');
if MainDir == 0
    return
end
animal = strfind(MainDir,'.');
animal = MainDir(animal-1:end);     % find animal number
root = dir(MainDir);                % list directory
slices = cell(0,1);
for n = 1:numel(root)
    sl = textscan(root(n).name ,'Slice %s');    %search for Slice directories
    if ~isempty(sl{1})
        slices{end+1} = sl{1}{1};               %remember Slice number
    end
end

ic = {'ipsi','contra'};
for sl = 1:numel(slices)                        %go through slice directories
    
    tracing = dir(fullfile(MainDir,sprintf('Slice %s',slices{sl}),'Tracings')); % list tracing folder
    tracing = {tracing(~cat(1,tracing.isdir)).name};    %chooses only file names
    for ci = 1:2
        %         all_trees.(ic{ci}) = cell(0,1);
        if exist(fullfile(MainDir,sprintf('Slice %s',slices{sl}),ic{ci},'PartStart.txt'),'file')
            if strfind(options,'-s')
                f= figure('Name',sprintf('RV-GFP Rat %s - Slice %s %s',animal,slices{sl},ic{ci}),'NumberTitle','off');
                set(gca,'YDir','reverse')
                hold on
            end
            txt = fileread(fullfile(MainDir,sprintf('Slice %s',slices{sl}),ic{ci},'PartStart.txt'));
            txt = textscan(txt,'%s','Delimiter','\n');
            txt = txt{1};
            part = cell(numel(txt),3);   %initialize
            ct=0;
            for p = 1:numel(txt)
                part(p,:) = textscan(txt{p},'%s %f64 %f64');    % scans the pth line
                
                ind = find(~cellfun(@isempty,strfind(tracing,sprintf('%s_%s_%s_part%s',animal,slices{sl},ic{ci},part{p,1}{1}))));
                
                if ~isempty(ind)
                    tree = load_tree(fullfile(MainDir,sprintf('Slice %s',slices{sl}),'Tracings',tracing{ind}));
                    if ~isstruct(tree{1}) && iscell(tree{1})
                       tree = tree{1}; 
                    end
                    for t = 1:numel(tree)
                        if strcmp(tree{t}.done,'yes')
                            tr = tree{t};
                            tr.X = tr.X + part{p,2} * tr.x_scale;
                            tr.Y = tr.Y + part{p,3} * tr.y_scale;
                            if strfind(options,'-s')
                                plot_tree(tr,[],[],[],[],'-p');
                            end
                            all_trees(str2num(slices{sl})).(ic{ci}){ct+t} = tr;
                        else
                            ct = ct -1;
                        end
                    end
                    ct = ct + t;
                end
            end
            if isfield(all_trees(str2num(slices{sl})),ic{ci}) && ~isempty (all_trees(str2num(slices{sl})).(ic{ci}))
                save_tree(all_trees(str2num(slices{sl})).(ic{ci}),fullfile(MainDir,sprintf('Slice %s',slices{sl}),sprintf('%s_%s_%s_all_trees.mtr',animal,slices{sl},ic{ci})))
            elseif strfind(options,'-s')
                delete(f)
            end
        end
    end
end