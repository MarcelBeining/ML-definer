%%
function gui_getMLcontours_big(src, event, action,par)

global F M sizM Z lines trees status contours params

if nargin > 3
    params = par;
elseif isempty(params)
   params.MLmode = 1;
   params.intpoints = 2000;
   params.relsmooth = 0.02;
   params.Zstep = 5;
end
if nargin<1
    action = 'init';
end

trees = [];
col = {'w','g','b','y','r'};


switch action,
    case 'move'
        if status(Z)~=6
            update_plot(F,sizM,Z,lines,contours,status,1)
        end
    case 'addpoint'
        if status(Z)~=6
            if isfield(lines{Z,status(Z)},'Vertices') && ~isempty(lines{Z,status(Z)}.Vertices)
                Point = getPoint(sizM,0);
                lines{Z,status(Z)}.Vertices = [lines{Z,status(Z)}.Vertices; Point(1), Point(2)];
            else
                lines{Z,status(Z)}.Vertices = getPoint(sizM,0);
            end
            update_plot(F,sizM,Z,lines,contours,status,2)
        end
    case 'deletepoint'
        if status(Z) ~=6 && isfield(lines{Z,status(Z)},'Vertices') && ~isempty(lines{Z,status(Z)}.Vertices)
            Point = getPoint(sizM,0);
            dist = sqrt((lines{Z,status(Z)}.Vertices(:,1)-Point(1)).^2+(lines{Z,status(Z)}.Vertices(:,2)-Point(2)).^2);
            [m indy] = min(dist);
            lines{Z,status(Z)}.Vertices(indy,:) = [];
            update_plot(F,sizM,Z,lines,contours,status,2)
        end
    case 'keypress'
        ui = 0;
        
        if strcmp(event.Character,'f') && Z < size(M,3)
            Z = Z + 1;
            ui = 3;
        elseif strcmp(event.Character,'d') && Z > 1
            Z = Z - 1 ;
            ui = 3;
        elseif strcmp(event.Character,' ')
            if size(lines{Z,status(Z)}.Vertices,1)<2
                errordlg('The contour needs to consist of at least two points','WindowStyle','modal' );
                return
            end
            if status(Z)~=6
                [lines{Z,status(Z)}.Vertices(:,1),lines{Z,status(Z)}.Vertices(:,2)] = poly2cw(lines{Z,status(Z)}.Vertices(:,1),lines{Z,status(Z)}.Vertices(:,2));
                lines{Z,status(Z)}.Faces = 1 : size(lines{Z,status(Z)}.Vertices,1);
            end
            if status(Z) == 5
                w = waitbar(0,'Calculating layer borders, please wait...');
               SGCL = interp_border(lines{Z,1}.Vertices,params.intpoints,1,'smooth',params.relsmooth); %!%!
               waitbar(1/6,w)
               fissura = interp_border(lines{Z,5}.Vertices,params.intpoints,1,'smooth',params.relsmooth); %!%!
               waitbar(2/6,w)
               GCL = interp_border(lines{Z,2}.Vertices,params.intpoints,1,'smooth',params.relsmooth); %!%!
               waitbar(3/6,w)
               lines{Z,2}.Vertices = GCL;
               lines{Z,5}.Vertices = fissura;
               lines{Z,1}.Vertices = SGCL;
               if params.MLmode == 1
               z = 0;
                    for o = 1: size(SGCL,1)
                        pairs=[SGCL(o,:);fissura(o,:)];
                        
                        [x y] = polyxpoly(pairs(:,1),pairs(:,2),lines{Z,2}.Vertices(:,1),lines{Z,2}.Vertices(:,2));
                        if numel(x)<2 && ~isempty(x)
                            vec = diff(pairs,1,1);
                            vec = vec / norm(vec);
                            dis = norm([x-pairs(1,1) y-pairs(1,2)]);
                            lines{Z,3}.Vertices(o-z,:) = pairs(1,:) + vec*dis*2;
                            lines{Z,4}.Vertices(o-z,:) = lines{Z,3}.Vertices(o-z,:) + (pairs(2,:)-lines{Z,3}.Vertices(o-z,:))/2;
                            
                        else
                            z = z + 1;
                        end
                    end
               elseif params.MLmode == 2
                   lines{Z,3}.Vertices = GCL + (fissura - GCL) / 3;
                   lines{Z,4}.Vertices = GCL + (fissura - GCL)* 2 / 3;
               end
                    waitbar(5/6,w)
                lines{Z,3}.Faces = 1:size(lines{Z,3}.Vertices,1);
                lines{Z,4}.Faces = lines{Z,3}.Faces;
                contours(Z*params.Zstep,:) = makeContours(lines,Z);
                waitbar(6/6,w)
                close(w)
            end
            m = [2 5 0 0 6 6];
            status(Z) = m(status(Z));
            ui=2;
        elseif strcmp(event.Key,'backspace')
            m = [1 1 0 0 2 5];
            ui = 2;

            lines{Z,status(Z)} = [];
            if status(Z) == 6
                contours(Z,:)=cell(1,4);
                lines{Z,3} = [];
                lines{Z,4} = [];
%                 lines{Z,3}.Vertices([1 end],:) = [];
%                 lines{Z,3}.Faces = lines{Z,3}.Faces(1:end-2);
%                 lines{Z,4}.Vertices([1 end],:) = [];
%                 lines{Z,4}.Faces = lines{Z,4}.Faces(1:end-2);
            end
            status(Z) = m(status(Z));
        elseif strcmp(event.Character,'+')
            contrast = caxis;
            contrast(2) = max(contrast(2)-10,contrast(1)+1);
            caxis(contrast);
        elseif strcmp(event.Character,'-')
            contrast = caxis;
            contrast(2) = min(contrast(2)+10,255);
            caxis(contrast);
        elseif strcmp(event.Character,'l')
            if exist(fullfile(params.PathName,'ML-contours_interpolated.mat'),'file')
                load(fullfile(params.PathName,'ML-contours_interpolated.mat'))
            end
        end
        update_plot(F,sizM,Z,lines,contours,status,ui,M)
        
    case 'init'
        [params.FileName,params.PathName] = uigetfile({'*.tif','*.tif - Multi-frame TIF file'},'Please open file');
        if params.FileName == 0
            return
        end
        info = imfinfo(fullfile(params.PathName,params.FileName));
        M = zeros(info(1).Height,info(1).Width,ceil(numel(info)/params.Zstep),sprintf('uint%d',ceil(info(1).BitDepth/8)*8));
        w = waitbar(0,'Loading Tif file...');
        zstep = 1:params.Zstep:numel(info);
        for n = 1:numel(zstep)
            M(:,:,n) = imread(fullfile(params.PathName,params.FileName),zstep(n));
            waitbar(n/size(M,3),w)
        end
        close(w)
        sizM = cat(2,[1;1;1],[size(M,1),size(M,2),size(M,3)]');
        sizM = [sizM([2 1],1);sizM([2 1],2)];
        contours = cell(numel(info),4);
        lines = cell(numel(info),5);
        Z = 1;
        
        status = ones(numel(info),1);
        F = figure('Name','ML-Definer 2013','NumberTitle','off');
        axis image;
        hold on
        set(gca,'YDir','reverse')
        colormap gray;
        update_plot(F,sizM,Z,lines,contours,status,3,M);
        
        %         set(gca,'ydir','reverse');
        
        %         set(gca,'visible','off');
        set(F,'KeyPressFcn',{@gui_getMLcontours_big, 'keypress'});
        set(F,'WindowButtonMotionFcn',{@gui_getMLcontours_big, 'move'});
        set(F,'WindowButtonDownFcn',{@gui_getMLcontours_big, 'addpoint'});
        set(F,'WindowScrollWheelFcn',{@gui_getMLcontours_big, 'deletepoint'});
        set(F,'CloseRequestFcn',{@gui_getMLcontours_big, 'close'});
        set(F,'Interruptible','off')
        
    case 'close'
        answer = questdlg('Save contours?');
        if strcmp(answer,'Yes')
            save(fullfile(params.PathName,'ML-contours.mat'),'contours')
        elseif strcmp(answer,'No')
            clearvars contours lines status
        elseif strcmp(answer,'Cancel') || isempty(answer)
            return
        end
        if exist('F','var') && ishandle(F)
            delete(F)
        end
        answer = questdlg('Proceed with contour processing and tree-MLyzing?');
        if strcmp(answer,'Yes')
            contours = interpz_contours(contours,[],'-ex-w');
            save(fullfile(params.PathName,'ML-contours_interpolated.mat'),'contours')
            answer = questdlg('Trees already transformed into whole DG?');
            if strcmp(answer,'No')
                TreeConverter
            end
            [FileName,PathName] = uigetfile({'*.mtr','*.mtr - Reconstructed Trees (Treestoolbox)'},'Please open whole DG trees file corresponding to the just segmentated DG');
            DGtrees = load_tree(fullfile(PathName,FileName));
            DGtrees = MLyzer_tree(DGtrees,contours,'-w');
            answer = questdlg('Plot result?');
            if strcmp(answer,'Yes')
                figure('Name','Result MLyzed trees','NumberTitle','off');
                hold on
                for t = 1:numel(DGtrees)
                   plot_tree(DGtrees{t},DGtrees{t}.R) 
                end
                for i = 1:4
                    patch(contours{1,i}.Vertices(:,1)*DGtrees{t}.x_scale,contours{1,i}.Vertices(:,2)*DGtrees{t}.y_scale,col{i},'FaceAlpha',0.4);
                end
            end
            save_tree(DGtrees,fullfile(PathName,sprintf('%s_MLyzed.mtr',FileName(1:end-4))))
        end
end

%%
    function update_plot(F,sizM,Z,lines,contours,status,mode,M)
        HP = get(F,'UserData');
%         hold off
        if mode == 1 && numel(HP)>1
            HHP = HP(2);
            delete(HHP(ishandle(HHP) & HHP~=0 ))
            HP(2)=0;
        elseif mode == 2
            HHP=HP(2:end);
            delete(HHP(ishandle(HHP) & HHP~=0 ))
            HP(2:end)=0;
            
        elseif mode == 3
            delete(HP(ishandle(HP) & HP~=0 ))
            HP=[];
            HP(1) = imagesc(M(:,:,Z));
            set(gca,'color','black')
%             xlim([1-sizM(3)/3 sizM(3)*4/3])
%             ylim([1-sizM(4)/3 sizM(4)*4/3])
        else 
            return
        end
        c = minmaxPoint(lines(Z,:));
        xlim([min(sizM(1),c(1)),max(sizM(3),c(2))])
        ylim([min(sizM(2),c(3)),max(sizM(4),c(4))])
%         if any(~ishandle(HP))
%             HP(~ishandle(HP))=0;
%         end
%         hold on
        switch status(Z)
            case 1
                title(gca,sprintf('Draw the SGZ/PML border line (plane %d)',Z*params.Zstep))
            case 2
                title(gca,sprintf('Now draw the GCL/IML border line (plane %d)',Z*params.Zstep))
                
            case 5
                title(gca,sprintf('Now draw the hippocampal fissure (plane %d)',Z*params.Zstep))
                
            case 6
                title(gca,sprintf('You are done with this plane (plane %d)',Z*params.Zstep))
        end
        if status(Z) == 6
            for i = 1:4
                HP(i+7)=patch(contours{Z*params.Zstep,i}.Vertices(:,1),contours{Z*params.Zstep,i}.Vertices(:,2),col{i});%,'FaceAlpha',0.3);%,'NextPlot', 'add');
%                 set(HP(i+1),'HitTest','off')
            end
        else
            if isfield(lines{Z,status(Z)},'Vertices') && ~isempty(lines{Z,status(Z)}.Vertices)
                Point = getPoint(sizM,0);
                if mode ~=1
                    HP(3) = plot(lines{Z,status(Z)}.Vertices(:,1),lines{Z,status(Z)}.Vertices(:,2),'r.');
                end
                HP(2) = plot([lines{Z,status(Z)}.Vertices(end,1);Point(1)],[lines{Z,status(Z)}.Vertices(end,2);Point(2)],sprintf('%s-',col{status(Z)}));
%             else
%                 Point = getPoint(sizM,0);
%                 HP(3) = plot(Point(1),Point(2),'r.');
            end
            if mode~=1
                for mo = 1:status(Z)
                    if isfield(lines{Z,mo},'Vertices') && ~isempty(lines{Z,mo}.Vertices)
                        HP(mo+3) = plot(lines{Z,mo}.Vertices(:,1),lines{Z,mo}.Vertices(:,2),sprintf('%s-',col{mo}));
                    end
                end
            end
        end
        set(F,'UserData',HP)
    end
%%

    function Point = getPoint(sizM,special)
        Point = get(gca,'CurrentPoint');
        Point = [Point(1),Point(3)];
%         Point = [max(sizM(1)-sizM(3)*4/3,Point(1)), max(sizM(1)-sizM(4)*4/3,Point(2))];
%         Point = [min(sizM(3)*4/3,Point(1)), min(sizM(4)*4/3,Point(2))];
        if special
            [nix ind]=min(abs(sizM-repmat(Point',[2 1])));
            Point(~rem(ind,2)+1)=sizM(ind);
        end
    end

%%
    function contour = makeContours(lines,Z)
        contour = cell(1,4);
        for i = 1:4
            line1 = lines{Z,i};
            line2 = lines{Z,i+1};

            line2.Vertices = flipud(line2.Vertices);

            contour{i}.Vertices = cat(1,line1.Vertices,line2.Vertices);%,line1.Vertices(1,:));
%             [contour{i}.Vertices(:,1),contour{i}.Vertices(:,2)] = poly2cw(contour{i}.Vertices(:,1),contour{i}.Vertices(:,2));
            contour{i}.Faces = 1:size(contour{i}.Vertices,1);
        end
    end
    function c = minmaxPoint(lins)
       nempty = find(~cellfun(@isempty,lins));
       c = [Inf Inf;-Inf -Inf];
       for ne = nempty
       c(1,:) = min([c(1,:);lins{ne}.Vertices]);
       c(2,:) = max([c(2,:);lins{ne}.Vertices]);
       end
    end
%%
end