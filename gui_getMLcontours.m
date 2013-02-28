function gui_getMLcontours(src, event, action)

global F M sizM Z lines trees mode contours

if nargin<1,
    action = 'init';
end

trees = [];
col = {'w','g','b','y','r'};


switch action,
    case 'move'
        if mode(Z)~=6
            update_plot(F,sizM,Z,lines,contours,mode)
        end
    case 'addpoint'
        if mode(Z)~=6
            if isfield(lines{Z,mode(Z)},'Vertices') && ~isempty(lines{Z,mode(Z)}.Vertices)
                Point = getPoint(sizM,0);
                lines{Z,mode(Z)}.Vertices = [lines{Z,mode(Z)}.Vertices; Point(1), Point(2)];
            else
                lines{Z,mode(Z)}.Vertices = getPoint(sizM,1);
            end
            update_plot(F,sizM,Z,lines,contours,mode)
        end
    case 'deletepoint'
        if mode(Z) ~=6 && isfield(lines{Z,mode(Z)},'Vertices') && ~isempty(lines{Z,mode(Z)}.Vertices)
            Point = getPoint(sizM,0);
            dist = sqrt((lines{Z,mode(Z)}.Vertices(:,1)-Point(1)).^2+(lines{Z,mode(Z)}.Vertices(:,2)-Point(2)).^2);
            [m indy] = min(dist);
            lines{Z,mode(Z)}.Vertices(indy,:) = [];
            update_plot(F,sizM,Z,lines,contours,mode)
        end
    case 'keypress'
        ui = 0;
        
        if strcmp(event.Character,'f') && Z < size(M,3)
            Z = Z + 1;
            ui = 2;
        elseif strcmp(event.Character,'d') && Z > 1
            Z = Z - 1 ;
            ui = 2;
        elseif strcmp(event.Character,' ')
            if ~any(lines{Z,mode(Z)}.Vertices(end,:) == 1) && ~any(lines{Z,mode(Z)}.Vertices(end,:) == sizM(3:4)')
                errordlg('The contour needs to end on the window border at both ends','WindowStyle','modal' );
                return
            elseif size(lines{Z,mode(Z)}.Vertices,1)<2
                errordlg('The contour needs to consist of at least two points','WindowStyle','modal' );
                return
            end
            if mode(Z)~=6
                lines{Z,mode(Z)}.Faces = 1 : size(lines{Z,mode(Z)}.Vertices,1);
            end
            if mode(Z) == 5
                th = linspace(0,2*pi,500)';

                [c(3,1) c(3,2)] = circfit(lines{Z,5}.Vertices(:,1),lines{Z,5}.Vertices(:,2));
                [c(2,1) c(2,2)] = circfit(lines{Z,2}.Vertices(:,1),lines{Z,2}.Vertices(:,2));
                [c(1,1) c(1,2)] = circfit(lines{Z,1}.Vertices(:,1),lines{Z,1}.Vertices(:,2));
                c = mean(c,1);
                test = [mean(lines{Z,5}.Vertices,1); c];
                if ~isempty(polyxpoly(test(:,1),test(:,2),lines{Z,2}.Vertices(:,1),lines{Z,2}.Vertices(:,2)))
                    r(1) = abs(mean(sqrt(sum((lines{Z,1}.Vertices - repmat(c,[size(lines{Z,1}.Vertices,1) 1])).^2,2))));
                    r(2) = abs(mean(sqrt(sum((lines{Z,2}.Vertices - repmat(c,[size(lines{Z,2}.Vertices,1) 1])).^2,2))));
                    r(5) = abs(mean(sqrt(sum((lines{Z,5}.Vertices - repmat(c,[size(lines{Z,5}.Vertices,1) 1])).^2,2))));
                    r(3) = abs((r(2)-r(1))*2 + r(1));
                    r(4) = abs((r(5)-r(3))/2 + r(3));
                    for ward = 1:2
                        lines{Z,ward+2}.Vertices = [c(1) + r(ward+2)*cos(th) c(2) + r(ward+2)*sin(th)];
                        lines{Z,ward+2}.Vertices = lines{Z,ward+2}.Vertices(lines{Z,ward+2}.Vertices(:,1) > sizM(1) & lines{Z,ward+2}.Vertices(:,2) >sizM(2) & lines{Z,ward+2}.Vertices(:,1) < sizM(3) & lines{Z,ward+2}.Vertices(:,2) < sizM(4),:);
                    end
                else
                    fissura = interpxy_contour(lines{Z,5}.Vertices,40); %!%!
                    SGCL = interpxy_contour(lines{Z,1}.Vertices,40); %!%!
                    
                    for o = 1: size(SGCL,1)
                        [nix,ind] = min(sqrt(sum((fissura-repmat(SGCL(o,:),[size(fissura,1) 1])).^2,2)));
                        pairs=[SGCL(o,:);fissura(ind,:)];
                        [x y] = polyxpoly(pairs(:,1),pairs(:,2),lines{Z,2}.Vertices(:,1),lines{Z,2}.Vertices(:,2));
                        if numel(x)<2
                            vec = diff(pairs,1,1);
                            vec = vec / norm(vec);
                            dis = norm([x-pairs(1,1) y-pairs(1,2)]);
                            lines{Z,3}.Vertices(o,:) = pairs(1,:) + vec*dis*2;
                            lines{Z,4}.Vertices(o,:) = lines{Z,3}.Vertices(o,:) + (pairs(2,:)-lines{Z,3}.Vertices(o,:))/2;
                            
                        end
                        
                    end
                    
                    
                    
                end
                %                 [x2 y2] = circfit(lines{Z,5}.Vertices(:,1),lines{Z,5}.Vertices(:,2));
%                 for o = 1: size(lines{Z,5}.Vertices,1)  %old algorithmus, find with middlepoint of circle
%                     pairs = [xc yc;lines{Z,5}.Vertices(o,: )];
%                     [x y] = polyxpoly(pairs(:,1),pairs(:,2),lines{Z,1}.Vertices(:,1),lines{Z,1}.Vertices(:,2));
%                     if numel(x)<2 && ~isempty(x)
%                         pairs(1,:) = [x y];
%                     else
%                         continue
%                     end
%                     [x y] = polyxpoly(pairs(:,1),pairs(:,2),lines{Z,2}.Vertices(:,1),lines{Z,2}.Vertices(:,2));
%                     if numel(x)<2 && ~isempty(x)
%                         vec = diff(pairs,1,1);
%                         vec = vec / norm(vec);
%                         dis = norm([x-pairs(1,1) y-pairs(1,2)]);
%                         lines{Z,3}.Vertices(o,:) = pairs(1,:) + vec*dis*2;
%                         lines{Z,4}.Vertices(o,:) = lines{Z,3}.Vertices(o,:) + (pairs(2,:)-lines{Z,3}.Vertices(o,:))/2;
%                         
%                     end
%                 end
                box = [1 1;sizM(3) 1;sizM(3:4)';1 sizM(4);1 1];
                 
                for o = 3:4
                    vec = diff(lines{Z,o}.Vertices([2 1],:),1,1);
                    vec = vec/norm(vec);
                    pairs = [lines{Z,o}.Vertices(1,:);lines{Z,o}.Vertices(1,:)+vec*norm(sizM)];
                    [x y] = polyxpoly(pairs(:,1),pairs(:,2),box(:,1),box(:,2));
                    if ~all(lines{Z,o}.Vertices(1,:) == [x y])
                        lines{Z,o}.Vertices = cat(1,[x y],lines{Z,o}.Vertices);
                    end
                    vec = diff(lines{Z,o}.Vertices([end-1 end],:),1,1);
                    vec = vec/norm(vec);
                    pairs = [lines{Z,o}.Vertices(end,:);lines{Z,o}.Vertices(end,:)+vec*norm(sizM)];
                    [x y] = polyxpoly(pairs(:,1),pairs(:,2),box(:,1),box(:,2));
                    if ~isempty(x) && ~all(lines{Z,o}.Vertices(end,:) == [x y])
                        lines{Z,o}.Vertices = cat(1,lines{Z,o}.Vertices,[x y]);
                    end
                end
                lines{Z,3}.Faces = 1:size(lines{Z,3}.Vertices,1);
                lines{Z,4}.Faces = lines{Z,3}.Faces;
                contours(Z,:) = makeContours(lines,Z,box);
            end
            m = [2 5 0 0 6 6];
            mode(Z) = m(mode(Z));
            ui=1;
        elseif strcmp(event.Key,'backspace')
            m = [1 1 0 0 2 5];
            ui = 1;

            lines{Z,mode(Z)} = [];
            if mode(Z) == 6
                contours(Z,:)=cell(1,4);
                lines{Z,3} = [];
                lines{Z,4} = [];
%                 lines{Z,3}.Vertices([1 end],:) = [];
%                 lines{Z,3}.Faces = lines{Z,3}.Faces(1:end-2);
%                 lines{Z,4}.Vertices([1 end],:) = [];
%                 lines{Z,4}.Faces = lines{Z,4}.Faces(1:end-2);
            end
            mode(Z) = m(mode(Z));
        end
        if ui==2
            update_plot(F,sizM,Z,lines,contours,mode,M)
        elseif ui==1
            update_plot(F,sizM,Z,lines,contours,mode)
        end
    case 'init'
        [FileName,PathName] = uigetfile({'*.tif','*.tif - Multi-frame TIF file'},'Please open file');
        if FileName == 0
            return
        end
        info = imfinfo(fullfile(PathName,FileName));
        M = zeros(info(1).Height,info(1).Width,numel(info),sprintf('uint%d',ceil(info(1).BitDepth/8)*8));
        w = waitbar(0,'Loading Tif file...');
        for n = 1:numel(info)
            M(:,:,n) = imread(fullfile(PathName,FileName),n);
            waitbar(n/numel(info),w)
        end
        close(w)
        sizM = cat(2,[1;1;1],[size(M,1),size(M,2),size(M,3)]');
        sizM = [sizM([2 1],1);sizM([2 1],2)];
        contours = cell(numel(info),4);
        lines = cell(numel(info),5);
        Z = 1;
        
        mode = ones(numel(info),1);
        F = figure;
        axis image;
%         hold on;

        colormap gray;
        update_plot(F,sizM,Z,lines,contours,mode,M);
        
        %         set(gca,'ydir','reverse');
        
        %         set(gca,'visible','off');
        set(F,'KeyPressFcn',{@gui_getMLcontours, 'keypress'});
        set(F,'WindowButtonMotionFcn',{@gui_getMLcontours, 'move'});
        set(F,'WindowButtonDownFcn',{@gui_getMLcontours, 'addpoint'});
        set(F,'WindowScrollWheelFcn',{@gui_getMLcontours, 'deletepoint'});
        set(F,'CloseRequestFcn',{@gui_getMLcontours, 'close'});
        set(F,'Interruptible','off')
        
    case 'close'
        answer = questdlg('Save contours?');
        if strcmp(answer,'Yes')
            save('contours.mat','contours')
        elseif strcmp(answer,'No')
            clearvars contours lines mode
        elseif strcmp(answer,'Cancel') || isempty(answer)
            return
        end
        delete(F)
        
end

%%
    function update_plot(F,sizM,Z,lines,contours,mode,M)
        HP = get(F,'UserData');
        hold off
        if nargin > 6
            delete(HP(ishandle(HP) & HP~=0 ))
            HP(1) = imagesc(M(:,:,Z));
        else
            HHP=HP(2:end);
            delete(HHP(ishandle(HHP) & HHP~=0 ))
        end
        if ~isempty(HP) && numel(HP)~=1
            HP(2:end)=[];
        end
        hold on
        switch mode(Z)
            case 1
                title(gca,sprintf('Draw the SGZ/PML border line (plane %d)',Z))
            case 2
                title(gca,sprintf('Now draw the GCL/IML border line (plane %d)',Z))
                
            case 5
                title(gca,sprintf('Now draw the hippocampal fissure (plane %d)',Z))
                
            case 6
                title(gca,sprintf('You are done with this plane (plane %d)',Z))
        end
        if mode(Z) == 6
            for i = 1:4
                HP(i+1)=patch(contours{Z,i}.Vertices(:,1),contours{Z,i}.Vertices(:,2),col{i});%,'FaceAlpha',0.3);%,'NextPlot', 'add');
%                 set(HP(i+1),'HitTest','off')
            end
        else
            if isfield(lines{Z,mode(Z)},'Vertices') && ~isempty(lines{Z,mode(Z)}.Vertices)
                Point = getPoint(sizM,0);
                HP(2) = plot(lines{Z,mode(Z)}.Vertices(:,1),lines{Z,mode(Z)}.Vertices(:,2),'r.');
                HP(3) = plot([lines{Z,mode(Z)}.Vertices(:,1);Point(1)],[lines{Z,mode(Z)}.Vertices(:,2);Point(2)],sprintf('%s-',col{mode(Z)}));
            else
                Point = getPoint(sizM,1);
                HP(2) = plot(Point(1),Point(2),'r.');
            end
            for mo = 1:mode(Z)-1
                if isfield(lines{Z,mo},'Vertices') && ~isempty(lines{Z,mo}.Vertices)
                    HP(mo+3) = plot(lines{Z,mo}.Vertices(:,1),lines{Z,mo}.Vertices(:,2),sprintf('%s-',col{mo}));
                end
            end
        end
        set(F,'UserData',HP)
    end
%%

    function Point = getPoint(sizM,special)
        Point = get(gca,'CurrentPoint');
        Point = [Point(1),Point(3)];
        Point(Point<1)=1;
        Point = [min(sizM(3),Point(1)), min(sizM(4),Point(2))];
        if special
            [nix ind]=min(abs(sizM-repmat(Point',[2 1])));
            Point(~rem(ind,2)+1)=sizM(ind);
        end
    end

%%
    function contour = makeContours(lines,Z,box)
        contour = cell(1,4);
        for i = 1:4
            line1 = lines{Z,i};
            line2 = lines{Z,i+1};
            corners = [line1.Vertices(1,:);line1.Vertices(end,:);line2.Vertices(1,:);line2.Vertices(end,:)];
            
            
            
            corner = cell(2,1);
            flip = false;
            for u = 1:2
                if any(diff(corners([u 3],:),1,1)==0)
                    if u ==1
                        flip = true;
                    end
                    corner{u} = [];
                elseif any(diff(corners([u 4],:),1,1)==0)
                    if u ==2
                        flip = true;
                    end
                    corner{u} = [];
                else
                    [nix ind] = min([pdist(corners([u 3],:)),pdist(corners([u 4],:))]);
                    if ind-u == 0
                        flip = true;
                    end
                    if ~isempty(intersect([corners(u,1),corners(2+ind,2)],box,'rows'))
                        corner{u} = [corners(u,1),corners(2+ind,2)];
                    else
                        corner{u} = [corners(2+ind,1),corners(u,2)];
                    end
                end
            end
            if flip
%                 line2.Faces = fliplr(line2.Faces);
                line2.Vertices = flipud(line2.Vertices);
            end
            contour{i}.Vertices = cat(1,line1.Vertices,corner{2},line2.Vertices,corner{1});%,line1.Vertices(1,:));
            [contour{i}.Vertices(:,1),contour{i}.Vertices(:,2)] = poly2cw(contour{i}.Vertices(:,1),contour{i}.Vertices(:,2));
            contour{i}.Faces = 1:size(contour{i}.Vertices,1);
        end
    end
%%
end