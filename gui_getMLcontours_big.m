%%
function gui_getMLcontours_big(src, event, action,par)

global F M sizM Z lines lines_corrected trees status contours params        %initialization

if nargin > 3
    params = par;
elseif isempty(params)          % then standard parameters are loaded
    params.MLmode = 2;           % Modes: 1) IML=GCL, OML,MML = Rest/2; 2) Rest /3 = IML,MML,OML
    params.intpoints = 100;     % Amount of points on drawn lines after interpol
    params.intdistance = 25;
    params.relsmooth = 0.02;     % filter size of smoothing, 0.02 = 2% of all points
    params.Zstep = 5;            % only every xth zplane is loaded (for big stacks)
    params.sample_rate = 1;
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
                lines{Z,status(Z)}.sample_rate = params.sample_rate;
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
        
        if strcmp(event.Character,'p')
            gui_getMLcontours_big([],[],'proceed');            
        elseif strcmp(event.Character,'f') && Z < size(M,3)         % move one plane forward
            Z = Z + 1;
            ui = 3;
        elseif strcmp(event.Character,'d') && Z > 1             % move one plane backward
            Z = Z - 1 ;
            ui = 3;
        elseif strcmp(event.Character,' ')                      % finish current line draw
            if status(Z) ~= 6 && size(lines{Z,status(Z)}.Vertices,1)<2
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
                
                curv = sqrt(sum(diff([SGCL GCL fissura],2,1).^2,2));
                
                SGCL = interp_border(SGCL,params.intpoints,3,[],[],curv); % interpolating SGCL line
                waitbar(1/6,w)
                fissura = interp_border(fissura,params.intpoints,3,[],[],curv); % interpolating GCL line
                waitbar(2/6,w)
                GCL = interp_border(GCL,params.intpoints,3,[],[],curv); % interpolating fissura line
                waitbar(3/6,w)
                
                lines{Z,2}.Vertices = GCL;
                lines{Z,2}.Faces = 1:params.intpoints;
                lines{Z,5}.Vertices = fissura;
                lines{Z,5}.Faces = 1:params.intpoints;
                lines{Z,1}.Vertices = SGCL;
                lines{Z,1}.Faces = 1:params.intpoints;
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
%                 lines{Z,3}.Vertices = interp_border(lines{Z,3}.Vertices,params.intdistance,2); % reducing size of line
%                 lines{Z,4}.Vertices = interp_border(lines{Z,4}.Vertices,params.intdistance,2); % reducing size of line
                lines{Z,3}.Faces = 1:size(lines{Z,3}.Vertices,1);   % make the connection array
                lines{Z,3}.sample_rate = params.sample_rate;
                lines{Z,4}.Faces = lines{Z,3}.Faces;                % "
                lines{Z,4}.sample_rate = params.sample_rate;
%                 contours(Z*params.Zstep,:) = makeContours(lines,Z); % construct the layer contours from the border lines
                waitbar(6/6,w)
                close(w)
            end
            m = [2 5 0 0 6 6];             % the way in which the mode changes from one to next
            status(Z) = m(status(Z));
            ui=2;
        elseif strcmp(event.Key,'backspace')    % go back one mode
            m = [1 1 0 0 2 5];             % the way in which the mode changes backwards
            ui = 2;
            

            if status(Z) == 6              % finish mode
                contours(Z,:)=cell(1,5);    % delete the created contours
                lines{Z,3} = [];            % delete IML border
                lines{Z,4} = [];            % delete MML/OML border
            else
                lines{Z,status(Z)} = [];       % the border line of this mode is deleted
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
%                 contours(setdiff(1:size(contours,1),1:params.Zstep:size(contours,1)),:)=[];     %deletes all contours that cannot be seen anyway (due to Zstep)
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
        if numel(info) == 1
%             try
%                 str = regexp(info.ImageDescription,'\n','split');
%                 strind = ~cellfun(@(x) isempty(strfind(x,'slices')),str);
%                 n = textscan(str{strind},'slices=%d');
%                 info = repmat(info,n{1},1);
%             catch
                answer = questdlg('Image is no stack or unknown format. Continue?','Yes','No');
                if ~strcmp(answer,'Yes')
                    return
                end
%             end
        end
        answer = inputdlg('Please give the sampling rate of the stack','Sampling Rate',1,{'1'});
        sample_rate = str2num(answer{1});
        if ~isempty(sample_rate) && isnumeric(sample_rate)
            params.sample_rate = sample_rate;
        end
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
        contours = cell(numel(info),5);         % initialize contours
        lines = cell(size(M,3),5);            % initialize border lines
        lines_corrected = cell(size(contours,1),5);
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
            save(fullfile(params.PathName,'ML-contours.mat'),'contours','lines','status')
%         elseif strcmp(answer,'No')
%             clearvars contours lines status
        elseif strcmp(answer,'Cancel') || isempty(answer)
            return
        end
        if exist('F','var') && ishandle(F)      % close figure
            delete(F)
        end
        
    case 'proceed'
        answer = questdlg('Proceed with contour processing and tree-MLyzing?');
        if strcmp(answer,'Yes')
            %             contours = interpz_contours(contours,[],'-ex-w');   % inter- and extrapolate contours to all zplanes
            if ~isempty(lines_corrected)
                answer = questdlg('Interpolate z-planes again?');
                if strcmp(answer,'Yes')
                    lines_corrected = cell(size(contours,1),5);
                    if size(lines,2) == 6
                        lines = lines(:,1:5);
                    end
                    nempty = find(all(~cellfun(@isempty,lines),2)); %find(status == 6);
                    %             nempty = find(~cellfun(@isempty,lines(:,5)));
                    lines_corrected(nempty * params.Zstep,:) = lines(nempty,:);
                    answer = questdlg('Extrapolate contours or use last available contour?','Extrapolation','Extrapolate','Last contour','Extrapolate');
                    if strcmp(answer,'Extrapolate')
                        lines_corrected = interpz_lines(lines_corrected,[],'-ex-w-s');   % inter- and extrapolate lines to all zplanes
                    else
                        lines_corrected = interpz_lines(lines_corrected,[],'-w-s');   % interpolate lines to all zplanes
                    end
                else
                    
                    col = {'k','y','b','g','r','w'};
                    figure;hold on
                    for i=1:size(lines_corrected,1)
                        for mc = 1:5
                            this_line = lines_corrected{i,mc};
                            this_line.Vertices(:,3)=i;
                            line(this_line.Vertices(:,1),this_line.Vertices(:,2),this_line.Vertices(:,3),'Color',col{mc});
                            %         this_patch = lines{i,mc};
                            %         this_patch.FaceAlpha = 0.2;
                            %         this_patch.FaceColor = col{mc};
                            %         this_patch.Vertices(:,3)=i;
                            %         patch(this_patch)
                        end
                    end
                end
                contours = makeContours(lines_corrected,1:size(contours,1)); % construct the layer contours from the border lines
                statusx = status;
                lines(1:end-1,:) = lines_corrected([1:(size(lines,1)-1)]*params.Zstep,:);
                lines(end,:) = lines_corrected(end,:);
                status(:) = 6;
                save(fullfile(params.PathName,'ML-contours_interpolated.mat'),'contours','lines','lines_corrected','status')
%                 status = statusx;
                answer = questdlg('Trees already transformed into whole DG?');
                if strcmp(answer,'No')
                    TreeConverter(params.PathName)
                elseif strcmp(answer,'Cancel') || strcmp(answer,'')
                    return
                end
                [FileName,PathName] = uigetfile({'*.mtr','*.mtr - Reconstructed Trees (Treestoolbox)'},'Please open whole DG trees file corresponding to the just segmentated DG',params.PathName);
                if isempty(FileName)
                    return
                end
                DGtrees = load_tree(fullfile(PathName,FileName));
                if ~isstruct(DGtrees{1}) && iscell(DGtrees{1})
                    DGtrees = DGtrees{1};
                end
                DGtrees = MLyzer_tree(DGtrees,contours,params.sample_rate,'-w');       % assign different regions to branches dependent on their location in the ML
                answer = questdlg('Plot result?');
                if strcmp(answer,'Yes')
                    figure('Name','Result MLyzed trees','NumberTitle','off');
                    hold on
                    for t = 1:numel(DGtrees)
                        plot_tree(DGtrees{t},DGtrees{t}.R)
                    end
                    for i = 1:5
                        patch(contours{1,i}.Vertices(:,1)*DGtrees{t}.x_scale*params.sample_rate , contours{1,i}.Vertices(:,2)*DGtrees{t}.y_scale*params.sample_rate , repmat(DGtrees{t}.z_scale,[numel(contours{1,i}.Vertices(:,1)) 1]),col{i},'FaceAlpha',0.4);
                        patch(contours{end,i}.Vertices(:,1)*DGtrees{t}.x_scale*params.sample_rate , contours{end,i}.Vertices(:,2)*DGtrees{t}.y_scale*params.sample_rate , repmat(size(contours,1)*DGtrees{t}.z_scale,[numel(contours{end,i}.Vertices(:,1)) 1]),col{i},'FaceAlpha',0.4);
                    end
                end
                save_tree({DGtrees},fullfile(PathName,sprintf('%s_MLyzed.mtr',FileName(1:end-4))))    %save MLyzed trees
            end
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
        if Z==0 %status(Z) == 6       % finished status, draw only patches
            for p = 1:4
                HP(p+7)=patch(contours{Z*params.Zstep,p}.Vertices(:,1),contours{Z*params.Zstep,p}.Vertices(:,2),col{p});%,'FaceAlpha',0.3);%,'NextPlot', 'add');
            end
        else
            if status(Z) ~= 6 && isfield(lines{Z,status(Z)},'Vertices') && ~isempty(lines{Z,status(Z)}.Vertices)
                Point = getPoint(sizM,0);
                if mode ~=1
                    HP(3) = plot(lines{Z,status(Z)}.Vertices(:,1),lines{Z,status(Z)}.Vertices(:,2),'r.');   % replot points
                end
                HP(2) = plot([lines{Z,status(Z)}.Vertices(end,1);Point(1)],[lines{Z,status(Z)}.Vertices(end,2);Point(2)],sprintf('%s-',col{status(Z)}));    % plot line to mouse
            end
            if mode~=1              % also plot all other lines if they exist
                if status(Z) == 6
                    mz = 5;
                else
                    mz = status(Z);
                end
                for mo = 1:mz
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
        answer = inputdlg('Give another point for SGCL or press Cancel','Defining SGCL zone..',1,{'[1,1]'});
        while 1
            
            if ~isempty(answer) && all(size(str2num(answer{1})) == [1,2])
                point = str2num(answer{1});
                p = plot3(point(1,1),point(1,2),1,'ro','MarkerSize',25);
                answerold = answer;
                answer = inputdlg('This point ok? Else give another point or press Cancel','Defining SGCL zone..',1,answer);
                if all(answerold{1}==answer{1})
                    break
                else
                    delete(p)
                end
            end
            if isempty(answer)
                break
            end
        end
        contour = cell(numel(Z),5);                        % initialize
        for zz = 1:numel(Z)
            thisZ = Z(zz);
            for l = 1:4
                line1 = lines{thisZ,l};                     % current inner line
                line2 = lines{thisZ,l+1};                   % current outer line

                line2.Vertices = flipud(line2.Vertices);    % change order of second line

                contour{zz,l+1}.Vertices = cat(1,line1.Vertices,line2.Vertices);%,line1.Vertices(1,:));     % put together both contours
                %             [contour{l}.Vertices(:,1),contour{l}.Vertices(:,2)] = poly2cw(contour{l}.Vertices(:,1),contour{l}.Vertices(:,2));
                contour{zz,l+1}.Faces = 1:size(contour{zz,l+1}.Vertices,1);   % create point order
                contour{zz,l+1}.sample_rate = params.sample_rate;
            end

            if ~isempty(answer) && all(size(str2num(answer{1})) == [1,2])
                contour{zz,1}.Vertices = [lines{thisZ,1}.Vertices ; str2num(answer{1})];       %SGCL
            else
                contour{zz,1}.Vertices = lines{thisZ,1}.Vertices;       %SGCL
            end
            contour{zz,1}.Faces = 1:size(contour{zz,1}.Vertices,1);   % create point order
            contour{zz,1}.sample_rate = params.sample_rate;
        end
    end

    function c = minmaxPoint(lins)  % find minimum and maximum points of all lines
        nempty = find(~cellfun(@isempty,lins)); % those nonempty line entries
        c = [Inf Inf;-Inf -Inf];    % initialize minmax array
        for ne = nempty
            c(1,:) = min([c(1,:);lins{ne}.Vertices],[],1);   % find minimum
            c(2,:) = max([c(2,:);lins{ne}.Vertices],[],1);   % find maximum
        end
    end
%%
end