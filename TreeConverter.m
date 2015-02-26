function [all_trees,MainDir] = TreeConverter(path,options)
if nargin < 2
    options = [];
end
if nargin < 1
    path = [];
end
MainDir = uigetdir(path,'Select Root Directory of Animal');
if MainDir == 0
    return
end

% animal = strfind(MainDir,'.');
% animal = MainDir(animal-1:end);     % find animal number
% if ~isempty(strfind(animal,' '))
%     animal = animal(1:strfind(animal,' ')-1);
% end
animal = regexp(MainDir,' ','split');
animal = animal{~cellfun(@(x) isempty(strfind(x,'.')),animal)};
% animal = str2double(regexp(animal,'\.','split'));
slices = cell(0,1);
if ~isempty(strfind(MainDir,'Slice'))
    ind = strfind(MainDir,'Slice');
    slices{end+1} = MainDir(ind+6:ind+8);
              %remember Slice number
    MainDir = MainDir(1:ind-1);

else
    root = dir(MainDir);                % list directory
    for n = 1:numel(root)
        sl = textscan(root(n).name ,'Slice %s');    %search for Slice directories
        if ~isempty(sl{1}) && root(n).isdir
            slices{end+1} = sl{1}{1};               %remember Slice number
        end
    end
end
icc =  {'Ipsilateral','Contralateral'};
ic = {'ipsi','contra'};
for sl = 1:numel(slices)                        %go through slice directories
    
    tracing = dir(fullfile(MainDir,sprintf('Slice %s',slices{sl}),'Tracings')); % list tracing folder
    tracing = {tracing(~cat(1,tracing.isdir)).name};    %chooses only file names
    tracing = tracing(cellfun(@(x) strcmp(x(end-3:end),'.mtr'),tracing));   %choose only mtrs
    for ci = 1:2
        %         all_trees.(ic{ci}) = cell(0,1);
        flag = false;
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
            flag = true;
            for p = 1:numel(txt)
                part(p,:) = textscan(txt{p},'%s %f64 %f64');    % scans the pth line
                if isempty(strfind(part{p,1}{1},'art'))
                    part{p,1}{1} = strcat('part',part{p,1}{1});
                end
            end

        elseif exist(fullfile(MainDir,sprintf('Slice %s',slices{sl}),ic{ci},sprintf('%s_partROIs.zip',icc{ci})),'file')
            [cvsROIs] = ReadImageJROI(fullfile(MainDir,sprintf('Slice %s',slices{sl}),ic{ci},sprintf('%s_partROIs.zip',icc{ci})));
%             [x, y] = cellfun(@(x) deal(x.vnRectBounds(1), x.vnRectBounds(2)),cvsROIs);
            part = cell(numel(cvsROIs),3);   %initialize            
            flag = true;
            for p = 1:numel(cvsROIs)
                part(p,1) = textscan(cvsROIs{p}.strName,'Part%s');
                if isempty(part{p,1})
                    part(p,1) = textscan(cvsROIs{p}.strName,'%s');
                    part{p,1}{1} =  part{p,1}{1}(cell2mat(strfind( part{p,1},'art'))-1:end);
                end
                if isempty(strfind(part{p,1}{1},'art'))
                    part{p,1}{1} = strcat('part',part{p,1}{1});
                end
                part(p,2:3) =  {cvsROIs{p}.vnRectBounds(2),cvsROIs{p}.vnRectBounds(1)}; % height is y!
            end
        end
        
        if flag
            ct=0;
            for p = 1:size(part,1)
                ind = find(~cellfun(@isempty,strfind(tracing,sprintf('%s_%s_%s_part%s',animal,slices{sl},ic{ci},part{p,1}{1}))) | ~cellfun(@isempty,strfind(tracing,sprintf('%s_%s_%s_%s',animal,slices{sl},ic{ci},part{p,1}{1}))) | ~cellfun(@isempty,strfind(tracing,sprintf('%s_part%s',icc{ci},part{p,1}{1}))) | ~cellfun(@isempty,strfind(tracing,sprintf('%s_%s',icc{ci},part{p,1}{1}))));
%                 if isempty(ind)
%                     ind = find(~cellfun(@isempty,strfind(tracing,sprintf('%s_%s_%s_part%s',animal,slices{sl},icc{ci},part{p,1}{1}))));
%                 end
                if ~isempty(ind)
                    if numel(ind) > 1
                        ind = ind(~cellfun(@isempty,strfind(tracing(ind),icc{ci})) | ~cellfun(@isempty,strfind(tracing(ind),ic{ci})));
                        if numel(ind) >1
                            warndlg(sprintf('There are more than one mtr files in Slice %s %s which have "%s" in their name. Now the first file is chosen. If wrong, please check and rename',slices{sl},icc{ci},part{p,1}{1}),'Multiple Files','replace')
                        end
                    end
                    tree = load_tree(fullfile(MainDir,sprintf('Slice %s',slices{sl}),'Tracings',tracing{ind(1)}));
                    if ~isstruct(tree{1}) && iscell(tree{1})
                       tree = tree{1}; 
                    end
                    for t = 1:numel(tree)
                        if isfield(tree{t},'done') && strcmp(tree{t}.done,'yes')
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
            if ~exist('all_trees','var')
                warndlg(sprintf('No trees found for RV-GFP Rat %s - Slice %s %s',animal,slices{sl},ic{ci}))
%                 delete(f)
            elseif numel(all_trees) >= str2num(slices{sl}) && isfield(all_trees(str2num(slices{sl})),ic{ci}) && ~isempty (all_trees(str2num(slices{sl})).(ic{ci}))
                save_tree({all_trees(str2num(slices{sl})).(ic{ci})},fullfile(MainDir,sprintf('Slice %s',slices{sl}),sprintf('%s_%s_%s_all_trees.mtr',animal,slices{sl},ic{ci})))
            end
        end
    end
end
oldd = pwd;
cd(MainDir)
while 1
   answer = questdlg('Split any file into supra/infra?','Split Trees','Yes','No','No');
   if strcmp(answer,'Yes')
      Split_Tree_File
   else
       cd(oldd)
       break
   end    
end