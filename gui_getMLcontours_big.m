%%
function gui_getMLcontours_big(src, event, action,par)

global F M sizM Z lines trees status contours params        %initialization

if nargin > 3
    params = par;
elseif isempty(params)          % then standard parameters are loaded
    params.MLmode = 2;           % Modes: 1) IML=GCL, OML,MML = Rest/2; 2) Rest /3 = IML,MML,OML
    params.intpoints = 2000;     % Amount of points on drawn lines after interpol
    params.intdistance = 25;
    params.relsmooth = 0.02;     % filter size of smoothing, 0.02 = 2% of all points
    params.Zstep = 5;            % only every xth zplane is loaded (for big stacks)
end
if nargin<1
    action = 'init';
end

trees = [];
col = {'w','g','b','y','r'};    % color order for contours (inner to outer layer)


switch action,
    case 'move'                                             % moving mouse
        if status(Z)~=6
            update_plot(F,sizM,Z,lines,contours,status,1)   % update only line to mouse
        end
    case 'addpoint'
        if status(Z)~=6
            if isfield(lines{Z,status(Z)},'Vertices') && ~isempty(lines{Z,status(Z)}.Vertices)
                Point = getPoint(sizM,0);
                lines{Z,status(Z)}.Vertices = [lines{Z,status(Z)}.Vertices; Point(1), Point(2)];    %add current point to line
            else
                lines{Z,status(Z)}.Vertices = getPoint(sizM,0);     % add first point of line
            end
            update_plot(F,sizM,Z,lines,contours,status,2)       %update everything except image
        end
    case 'deletepoint'
        if status(Z) ~=6 && isfield(lines{Z,status(Z)},'Vertices') && ~isempty(lines{Z,status(Z)}.Vertices)
            Point = getPoint(sizM,0);
            dist = sqrt((lines{Z,status(Z)}.Vertices(:,1)-Point(1)).^2+(lines{Z,status(Z)}.Vertices(:,2)-Point(2)).^2);
            [m indy] = min(dist);                               % find nearest point of the line to mouse cursor
            lines{Z,status(Z)}.Vertices(indy,:) = [];           % delete this point
            update_plot(F,sizM,Z,lines,contours,status,2)       %update everything except image
        end
    case 'keypress'
        ui = 0;
        
        if strcmp(event.Character,'f') && Z < size(M,3)         % move one plane forward
            Z = Z + 1;
            ui = 3;
        elseif strcmp(event.Character,'d') && Z > 1             % move one plane backward
            Z = Z - 1 ;
            ui = 3;
        elseif strcmp(event.Character,' ')                      % finish current line draw
            if size(lines{Z,status(Z)}.Vertices,1)<2
                errordlg('The contour needs to consist of at least two points','WindowStyle','modal' );
                return
            end
            if status(Z)~=6                                     % line is finished and sorted clockwise
                [lines{Z,status(Z)}.Vertices(:,1),lines{Z,status(Z)}.Vertices(:,2)] = poly2cw(lines{Z,status(Z)}.Vertices(:,1),lines{Z,status(Z)}.Vertices(:,2));
                lines{Z,status(Z)}.Faces = 1 : size(lines{Z,status(Z)}.Vertices,1);
            end
            if status(Z) == 5
                w = waitbar(0,'Calculating layer borders, please wait...');
                SGCL = interp_border(lines{Z,1}.Vertices,params.intpoints,1,'smooth',params.relsmooth); % interpolating SGCL line
                waitbar(1/6,w)
                fissura = interp_border(lines{Z,5}.Vertices,params.intpoints,1,'smooth',params.relsmooth); % interpolating GCL line
                waitbar(2/6,w)
                GCL = interp_border(lines{Z,2}.Vertices,params.intpoints,1,'smooth',params.relsmooth); % interpolating fissura line
                waitbar(3/6,w)
                %                lines{Z,2}.Vertices = GCL;
                %                lines{Z,5}.Vertices = fissura;
                %                lines{Z,1}.Vertices = SGCL;
                if params.MLmode == 1
                    z = 0;
                    for o = 1: size(SGCL,1)
                        pairs=[SGCL(o,:);fissura(o,:)];     % place together one point of SGCL and fissura
                        
                        [x y] = polyxpoly(pairs(:,1),pairs(:,2),lines{Z,2}.Vertices(:,1),lines{Z,2}.Vertices(:,2)); % find intersection between SGCL to fissura points and GCL line
                        if numel(x)<2 && ~isempty(x)        % something is found
                            
                            vec = diff(pairs,1,1);          % make unit vector from SGCL to fissura
                            vec = vec / norm(vec);
                            dis = norm([x-pairs(1,1) y-pairs(1,2)]);    % measure distance from SGCL to GCL intersection
                            lines{Z,3}.Vertices(o-z,:) = pairs(1,:) + vec*dis*2;    % point on IML border is same dis to GCL as SGCL point
                            lines{Z,4}.Vertices(o-z,:) = lines{Z,3}.Vertices(o-z,:) + (pairs(2,:)-lines{Z,3}.Vertices(o-z,:))/2;    % rest distance betw IML and fissura is split in two halfths
                        else
                            z = z + 1;
                        end
                    end
                elseif params.MLmode == 2
                    lines{Z,3}.Vertices = GCL + (fissura - GCL) / 3;    % IML is just one third of distance GCL/fissura
                    lines{Z,4}.Vertices = GCL + (fissura - GCL)* 2 / 3; % same for MML/OML
                end
                waitbar(5/6,w)
                lines{Z,3}.Vertices = interp_border(lines{Z,3}.Vertices,params.intdistance,2); % reducing size of line
                lines{Z,4}.Vertices = interp_border(lines{Z,4}.Vertices,params.intdistance,2); % reducing size of line
                lines{Z,3}.Faces = 1:size(lines{Z,3}.Vertices,1);   % make the connection array
                lines{Z,4}.Faces = lines{Z,3}.Faces;                % "
                contours(Z*params.Zstep,:) = makeContours(lines,Z); % construct the layer contours from the border lines
                waitbar(6/6,w)
                close(w)
            end
            m = [2 5 0 0 6 6];             % the way in which the mode changes from one to next
            status(Z) = m(status(Z));
            ui=2;
        elseif strcmp(event.Key,'backspace')    % go back one mode
            m = [1 1 0 0 2 5];             % the way in which the mode changes backwards
            ui = 2;
            
            lines{Z,status(Z)} = [];       % the border line of this mode is deleted
            if status(Z) == 6              % finish mode
                contours(Z,:)=cell(1,4);    % delete the created contours
                lines{Z,3} = [];            % delete IML border
                lines{Z,4} = [];            % delete MML/OML border
            end
            status(Z) = m(status(Z));
        elseif strcmp(event.Character,'+')  % increase contrast
            contrast = caxis;
            contrast(2) = max(contrast(2)-10,contrast(1)+1);
            caxis(contrast);
        elseif strcmp(event.Character,'-')  % decrease contrast
            contrast = caxis;
            contrast(2) = min(contrast(2)+10,255);
            caxis(contrast);
        elseif strcmp(event.Character,'l')  % load a previously saved
            ui = 2;
            answer = questdlg('Load normal or interpolated contours?','Loading...','Normal','Interpolated','Cancel','Interpolated');
            if strcmp(answer,'Interpolated')
                if exist(fullfile(params.PathName,'ML-contours_interpolated.mat'),'file')
                    load(fullfile(params.PathName,'ML-contours_interpolated.mat'))
                else
                    errordlg('No previously saved interpolated contour file found')
                end
                contours(setdiff(1:size(contours,1),1:params.Zstep:size(contours,1)),:)=[];     %deletes all contours that cannot be seen anyway (due to Zstep)
            elseif strcmp(answer,'Normal')
                if exist(fullfile(params.PathName,'ML-contours.mat'),'file')
                    load(fullfile(params.PathName,'ML-contours.mat'))
                else
                    errordlg('No previously saved interpolated contour file found')
                end
            else
                ui = 0;
            end
        elseif strcmp(event.Character,'s')
            save(fullfile(params.PathName,'ML-contours.mat'),'contours','lines','status')
            msgbox('The contours have been saved','Save')
        end
        update_plot(F,sizM,Z,lines,contours,status,ui,M)
        
    case 'init'             % initialization at the beginning of the program
        [params.FileName,params.PathName] = uigetfile({'*.tif','*.tif - Multi-frame TIF file'},'Please open file');
        if params.FileName == 0     % if cancel o.ä. was pressed
            return
        end
        info = imfinfo(fullfile(params.PathName,params.FileName));  % get picture information
        M = zeros(info(1).Height,info(1).Width,ceil(numel(info)/params.Zstep),sprintf('uint%d',ceil(info(1).BitDepth/8)*8));    % construct picture matrix (zstep,BitDepth important)
        w = waitbar(0,'Loading Tif file...');
        zstep = 1:params.Zstep:numel(info);         % define the z-planes that are loaded
        for n = 1:numel(zstep)
            M(:,:,n) = imread(fullfile(params.PathName,params.FileName),zstep(n)); % read corresponding tiff-layer
            waitbar(n/size(M,3),w)
        end
        close(w)
        sizM = cat(2,[1;1;1],[size(M,1),size(M,2),size(M,3)]'); % construct matrix size info
        sizM = [sizM([2 1],1);sizM([2 1],2)];
        contours = cell(numel(info),4);         % initialize contours
        lines = cell(size(M,3),5);            % initialize border lines
        Z = 5;                                  % initialize current Z plane
        status = ones(size(M,3),1);           % initialize Mode array
        F = figure('Name','ML-Definer 2013','NumberTitle','off');   % initialize figure
        axis image;
        hold on
        set(gca,'YDir','reverse')
        colormap gray;
        
        update_plot(F,sizM,Z,lines,contours,status,3,M);    % start update function first time
        
        % initialize figure callback functions
        set(F,'KeyPressFcn',{@gui_getMLcontours_big, 'keypress'});
        set(F,'WindowButtonMotionFcn',{@gui_getMLcontours_big, 'move'});
        set(F,'WindowButtonDownFcn',{@gui_getMLcontours_big, 'addpoint'});
        set(F,'WindowScrollWheelFcn',{@gui_getMLcontours_big, 'deletepoint'});
        set(F,'CloseRequestFcn',{@gui_getMLcontours_big, 'close'});
        set(F,'Interruptible','off')
        
    case 'close'            % closing the program
        answer = questdlg('Save contours?');
        if strcmp(answer,'Yes')
            save(fullfile(params.PathName,'ML-contours.mat'),'contours','lines')
        elseif strcmp(answer,'No')
            clearvars contours lines status
        elseif strcmp(answer,'Cancel') || isempty(answer)
            return
        end
        if exist('F','var') && ishandle(F)      % close figure
            delete(F)
        end
        answer = questdlg('Proceed with contour processing and tree-MLyzing?');
        if strcmp(answer,'Yes')
            contours = interpz_contours(contours,[],'-ex-w');   % inter- and extrapolate contours to all zplanes 
            status(:) = 6;
            save(fullfile(params.PathName,'ML-contours_interpolated.mat'),'contours','lines','status')
            answer = questdlg('Trees already transformed into whole DG?');
            if strcmp(answer,'No')
                TreeConverter
            elseif strcmp(answer,'Cancel') || strcmp(answer,'')
               return 
            end
            [FileName,PathName] = uigetfile({'*.mtr','*.mtr - Reconstructed Trees (Treestoolbox)'},'Please open whole DG trees file corresponding to the just segmentated DG');
            DGtrees = load_tree(fullfile(PathName,FileName));
            DGtrees = MLyzer_tree(DGtrees,contours,'-w');       % assign different regions to branches dependent on their location in the ML
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
            save_tree(DGtrees,fullfile(PathName,sprintf('%s_MLyzed.mtr',FileName(1:end-4))))    %save MLyzed trees
        end
end

%%
    function update_plot(F,sizM,Z,lines,contours,status,mode,M)     % updates figure
        HP = get(F,'UserData');     % get handles to figure children
        %         hold off
        if mode == 1 && numel(HP)>1 % update only line to cursor (moving)
            HHP = HP(2);
            delete(HHP(ishandle(HHP) & HHP~=0 ))
            HP(2)=0;
        elseif mode == 2        % update all lines
            HHP=HP(2:end);
            delete(HHP(ishandle(HHP) & HHP~=0 ))
            HP(2:end)=0;
            
        elseif mode == 3        % also update image
            delete(HP(ishandle(HP) & HP~=0 ))
            HP=[];
            HP(1) = imagesc(M(:,:,Z));
            set(gca,'color','black')
        else
            return
        end
        c = minmaxPoint(lines(Z,:));    % get minimum and maximum point of all lines
        xlim([min(sizM(1),c(1)),max(sizM(3),c(2))]) % define image size correspondingly
        ylim([min(sizM(2),c(3)),max(sizM(4),c(4))])
        
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
        if status(Z) == 6       % finished status, draw only patches
            for p = 1:4
                HP(p+7)=patch(contours{Z*params.Zstep,p}.Vertices(:,1),contours{Z*params.Zstep,p}.Vertices(:,2),col{p});%,'FaceAlpha',0.3);%,'NextPlot', 'add');
            end
        else
            if isfield(lines{Z,status(Z)},'Vertices') && ~isempty(lines{Z,status(Z)}.Vertices)
                Point = getPoint(sizM,0);
                if mode ~=1
                    HP(3) = plot(lines{Z,status(Z)}.Vertices(:,1),lines{Z,status(Z)}.Vertices(:,2),'r.');   % replot points
                end
                HP(2) = plot([lines{Z,status(Z)}.Vertices(end,1);Point(1)],[lines{Z,status(Z)}.Vertices(end,2);Point(2)],sprintf('%s-',col{status(Z)}));    % plot line to mouse
            end
            if mode~=1              % also plot all other lines if they exist
                for mo = 1:status(Z)
                    if isfield(lines{Z,mo},'Vertices') && ~isempty(lines{Z,mo}.Vertices)
                        HP(mo+3) = plot(lines{Z,mo}.Vertices(:,1),lines{Z,mo}.Vertices(:,2),sprintf('%s-',col{mo}));
                    end
                end
            end
        end
        set(F,'UserData',HP)            % store handles
    end
%%

    function Point = getPoint(sizM,special)         % gives current point
        Point = get(gca,'CurrentPoint');
        Point = [Point(1),Point(3)];
        %         Point = [max(sizM(1)-sizM(3)*4/3,Point(1)), max(sizM(1)-sizM(4)*4/3,Point(2))];
        %         Point = [min(sizM(3)*4/3,Point(1)), min(sizM(4)*4/3,Point(2))];
        if special                  % point has to lie in the picture
            [nix ind]=min(abs(sizM-repmat(Point',[2 1])));
            Point(~rem(ind,2)+1)=sizM(ind);
        end
    end

%%
    function contour = makeContours(lines,Z)        % make layer contours by using the border lines
        contour = cell(1,4);                        % initialize
        for l = 1:4
            line1 = lines{Z,l};                     % current inner line
            line2 = lines{Z,l+1};                   % current outer line
            
            line2.Vertices = flipud(line2.Vertices);    % change order of second line
            
            contour{l}.Vertices = cat(1,line1.Vertices,line2.Vertices);%,line1.Vertices(1,:));     % put together both contours
            %             [contour{l}.Vertices(:,1),contour{l}.Vertices(:,2)] = poly2cw(contour{l}.Vertices(:,1),contour{l}.Vertices(:,2));
            contour{l}.Faces = 1:size(contour{l}.Vertices,1);   % create point order
        end
    end
    function c = minmaxPoint(lins)  % find minimum and maximum points of all lines
        nempty = find(~cellfun(@isempty,lins)); % those nonempty line entries
        c = [Inf Inf;-Inf -Inf];    % initialize minmax array
        for ne = nempty
            c(1,:) = min([c(1,:);lins{ne}.Vertices]);   % find minimum
            c(2,:) = max([c(2,:);lins{ne}.Vertices]);   % find maximum
        end
    end
%%
end