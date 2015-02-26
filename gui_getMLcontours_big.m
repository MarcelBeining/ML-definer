%%
function gui_getMLcontours_big(src, event, action,par)

global F M sizM Z lines lines_corrected trees status contours params endstatus       %initialization

if nargin > 3
    params = par;
elseif isempty(params)          % then standard parameters are loaded
    params.MLmode = 2;           % Modes: 1) IML=GCL, OML,MML = Rest/2; 2) Rest /3 = IML,MML,OML
    params.intpoints = 201;     % Amount of points on drawn lines after interpol
    params.intdistance = 25;
    params.relsmooth = 2 / params.intpoints ;     % filter size of smoothing, 0.02 = 2% of all points
    params.Zstep = 5;            % only every xth zplane is loaded (for big stacks)
    params.sample_rate = 1;
    params.curvfilt = 0.05;
    params.min_intervall = 201;     % min intervall when using curve interpolation with pointdistribution
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
%                 sempty = find(~cellfun(@isempty,lines(Z,:)));
%                 [m, indi] = min(cellfun(@(x) min(sqrt((x.Vertices(:,1)-Point(1)).^2+(x.Vertices(:,2)-Point(2)).^2)),lines(Z,sempty)));
%                 if size(lines{Z,sempty(indi)}.Vertices,1) == 1
%                     lines{Z,sempty(indi)}.Vertices = [lines{Z,sempty(indi)}.Vertices; Point(1), Point(2)];    %add current point to line
%                 else
%                     dist = sqrt((lines{Z,sempty(indi)}.Vertices(:,1)-Point(1)).^2+(lines{Z,sempty(indi)}.Vertices(:,2)-Point(2)).^2);
%                     [m, indy] = min(dist);    % find nearest point of the line to mouse cursor
%                     %                     dist(indy) = Inf;
%                     if indy+1 > numel(dist)
%                          lines{Z,sempty(indi)}.Vertices = [lines{Z,sempty(indi)}.Vertices; Point(1), Point(2)];
%                     elseif indy-1 == 0
%                         lines{Z,sempty(indi)}.Vertices = [ Point(1), Point(2); lines{Z,sempty(indi)}.Vertices];
%                     else
%                         [m, indy2] = min(dist([indy-1,indy+1]));
%                         if indy2 == 1
%                             lines{Z,sempty(indi)}.Vertices = [lines{Z,sempty(indi)}.Vertices(1:indy-1,:); Point(1), Point(2); lines{Z,sempty(indi)}.Vertices(indy:end,:)];    %add current point to line
%                         else
%                             lines{Z,sempty(indi)}.Vertices = [lines{Z,sempty(indi)}.Vertices(1:indy,:); Point(1), Point(2); lines{Z,sempty(indi)}.Vertices(indy+1:end,:)];    %add current point to line
%                         end
%                     end
%                 end
            else
                lines{Z,status(Z)}.Vertices = getPoint(sizM,0);     % add first point of line
                lines{Z,status(Z)}.sample_rate = params.sample_rate;
            end
            update_plot(F,sizM,Z,lines,contours,status,2)       %update everything except image
        end
    case 'deletepoint'
        if status(Z) ~=6 && isfield(lines{Z,status(Z)},'Vertices') && ~isempty(lines{Z,status(Z)}.Vertices)
            Point = getPoint(sizM,0);
%             dist = sqrt((lines{Z,status(Z)}.Vertices(:,1)-Point(1)).^2+(lines{Z,status(Z)}.Vertices(:,2)-Point(2)).^2);
            sempty = find(~cellfun(@isempty,lines(Z,:)));
            [m, indi] = min(cellfun(@(x) min(sqrt((x.Vertices(:,1)-Point(1)).^2+(x.Vertices(:,2)-Point(2)).^2)),lines(Z,sempty)));
            dist = sqrt((lines{Z,sempty(indi)}.Vertices(:,1)-Point(1)).^2+(lines{Z,sempty(indi)}.Vertices(:,2)-Point(2)).^2);
            [m, indy] = min(dist);    % find nearest point of the line to mouse cursor
            lines{Z,sempty(indi)}.Vertices(indy,:) = [];           % delete this point
            update_plot(F,sizM,Z,lines,contours,status,2)       %update everything except image
        end
    case 'keypress'
        ui = 0;
        
        if strcmp(event.Character,'p')
            gui_getMLcontours_big([],[],'proceed');
        elseif strcmp(event.Character,'e')
            if status(Z)~=6
                if isfield(lines{Z,status(Z)},'Vertices') && ~isempty(lines{Z,status(Z)}.Vertices)
                    Point = getPoint(sizM,0);
                    %                 lines{Z,status(Z)}.Vertices = [lines{Z,status(Z)}.Vertices; Point(1), Point(2)];    %add current point to line
                    sempty = find(~cellfun(@isempty,lines(Z,:)));
                    [m, indi] = min(cellfun(@(x) min(sqrt((x.Vertices(:,1)-Point(1)).^2+(x.Vertices(:,2)-Point(2)).^2)),lines(Z,sempty)));
                    if size(lines{Z,sempty(indi)}.Vertices,1) == 1
                        lines{Z,sempty(indi)}.Vertices = [lines{Z,sempty(indi)}.Vertices; Point(1), Point(2)];    %add current point to line
                    else
                        dist = sqrt((lines{Z,sempty(indi)}.Vertices(:,1)-Point(1)).^2+(lines{Z,sempty(indi)}.Vertices(:,2)-Point(2)).^2);
                        [m, indy] = min(dist);    % find nearest point of the line to mouse cursor
                        %                     dist(indy) = Inf;
                        if indy+1 > numel(dist)
                            lines{Z,sempty(indi)}.Vertices = [lines{Z,sempty(indi)}.Vertices; Point(1), Point(2)];
                        elseif indy-1 == 0
                            lines{Z,sempty(indi)}.Vertices = [ Point(1), Point(2); lines{Z,sempty(indi)}.Vertices];
                        else
                            [m, indy2] = min(dist([indy-1,indy+1]));
                            if indy2 == 1
                                lines{Z,sempty(indi)}.Vertices = [lines{Z,sempty(indi)}.Vertices(1:indy-1,:); Point(1), Point(2); lines{Z,sempty(indi)}.Vertices(indy:end,:)];    %add current point to line
                            else
                                lines{Z,sempty(indi)}.Vertices = [lines{Z,sempty(indi)}.Vertices(1:indy,:); Point(1), Point(2); lines{Z,sempty(indi)}.Vertices(indy+1:end,:)];    %add current point to line
                            end
                        end
                    end
                else
                    lines{Z,status(Z)}.Vertices = getPoint(sizM,0);     % add first point of line
                    lines{Z,status(Z)}.sample_rate = params.sample_rate;
                end
                update_plot(F,sizM,Z,lines,contours,status,2)       %update everything except image
            end
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
%                 [lines{Z,status(Z)}.Vertices(:,1),lines{Z,status(Z)}.Vertices(:,2)] = poly2cw(lines{Z,status(Z)}.Vertices(:,1),lines{Z,status(Z)}.Vertices(:,2));
                lines{Z,status(Z)}.Faces = 1 : size(lines{Z,status(Z)}.Vertices,1);
            end
            if status(Z) == 5
                w = waitbar(0,'Calculating layer borders, please wait...');
                blade = {'','_supra','_infra'};
                save(fullfile(params.PathName,sprintf('ML-contours%s.mat',blade{params.blade+1})),'contours','lines','status')
             
                SGCL = interp_border(lines{Z,1}.Vertices,params.intpoints,1,'smooth',params.relsmooth); % interpolating SGCL line
%                 SGCL = lines{Z,1}.Vertices;
                waitbar(1/6,w)
                fissura = interp_border(lines{Z,5}.Vertices,params.intpoints,1,'smooth',params.relsmooth); % interpolating GCL line
                waitbar(2/6,w)
                GCL = interp_border(lines{Z,2}.Vertices,params.intpoints,1,'smooth',params.relsmooth); % interpolating fissura line
                waitbar(3/6,w)
                
%                %%old  curv = [0; sqrt(sum(diff([SGCL GCL fissura],2,1).^2,2))];    % make pdistr for interp dependent on curvature level
%                 curv = [0;sqrt(sum(diff(SGCL,2,1).^2,2));0];    % make pdistr for interp dependent on curvature level
%                 curv = [0;sqrt(sum(diff(GCL,2,1).^2,2));0];    % make pdistr for interp dependent on curvature level
                curv = sqrt(sum(diff(SGCL,2,1).^2,2));    % make pdistr for interp dependent on curvature level
%                 curv=curv.^2; % enhance biggest curves
                if round(size(curv,1)*params.curvfilt) ~= size(curv,1)*params.curvfilt
                    kern = rem(floor(size(curv,1)*params.curvfilt),2) * floor(size(curv,1)*params.curvfilt) + rem(ceil(size(curv,1)*params.curvfilt),2) * ceil(size(curv,1)*params.curvfilt); % rounds to the odd number
                else
                    kern = size(curv,1)*params.curvfilt + ~rem(size(curv,1)*params.curvfilt ,2) ;
                end
                tmp = curv([1,end]);
                curv = convn (padarray(curv,(kern-1)/2,'replicate'), ones (1, kern)' / kern, 'valid');
                curv = flipud(convn (padarray(flipud(curv),(kern-1)/2,'replicate'), ones (1, kern)' / kern, 'valid'));
                curv([1,end]) = tmp;   % accounts for smoothing error at edges
                
                [nix, ind] = max(curv);
                curv = (1:0.1:10).^2 ;  % earlier exp
                tSGCL = SGCL;
%                 tSGCL(1:ceil(params.intpoints/2),:) = interp_border(SGCL(1:ind+1,:),ceil(params.intpoints/2),3,[],[],curv(1:ind),params.min_intervall); % interpolating SGCL line
%                 tSGCL(ceil(params.intpoints/2)+1:end,:) = interp_border(SGCL(ind+2:end,:),floor(params.intpoints/2),3,[],[],curv(ind+1:end),params.min_intervall); % interpolating SGCL line
                tSGCL(1:ceil(params.intpoints/2),:) = interp_border(SGCL(1:ind+1,:),ceil(params.intpoints/2),3,[],[],curv,params.min_intervall); % interpolating SGCL line
                curv = fliplr(curv);
%                 curv(1) = 0;
                tSGCL(ceil(params.intpoints/2)+1:end,:) = interp_border(SGCL(ind+2:end,:),floor(params.intpoints/2),4,[],[],curv,params.min_intervall); % interpolating SGCL line
                SGCL = tSGCL;
%                 GCL = interp_border(GCL,params.intpoints,3,[],[],curv,params.min_intervall); % interpolating SGCL line
%                 waitbar(1/6,w)
                
                [GCL fissura] = make_orthogonal_borders(SGCL,GCL,fissura);
%                 [SGCL fissura] = make_orthogonal_borders(GCL,SGCL,fissura);     %CAUTION Changed to gcl
%             %%old    fissura = interp_border(fissura,params.intpoints,3,[],[],curv); % interpolating GCL line
%             %%old    waitbar(2/6,w)
%             %%old    GCL = interp_border(GCL,params.intpoints,3,[],[],curv); % interpolating fissura line
%             %%old    waitbar(3/6,w)
                
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
            endstatus = false;
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
            blade = {'','_supra','_infra'};
            
            answer = questdlg('Load normal or interpolated contours?','Loading...','Normal','Interpolated','Cancel','Interpolated');
            if strcmp(answer,'Interpolated')
                endstatus = true;
                if exist(fullfile(params.PathName,sprintf('ML-contours_interpolated%s.mat',blade{params.blade+1})),'file')
                    load(fullfile(params.PathName,sprintf('ML-contours_interpolated%s.mat',blade{params.blade+1})))
                else
                    errordlg('No previously saved interpolated contour file for this stack found')
                end
%                 contours(setdiff(1:size(contours,1),1:params.Zstep:size(contours,1)),:)=[];     %deletes all contours that cannot be seen anyway (due to Zstep)
            elseif strcmp(answer,'Normal')
                endstatus = false;
                if exist(fullfile(params.PathName,sprintf('ML-contours%s.mat',blade{params.blade+1})),'file')
                    load(fullfile(params.PathName,sprintf('ML-contours%s.mat',blade{params.blade+1})))
                    status = sum(~cellfun(@isempty,lines),2) +1;
                    status(status==4) = 5;
                else
                    errordlg('No previously saved interpolated contour file for this stack found')
                end
            else
                ui = 0;
            end
        elseif strcmp(event.Character,'s')
            blade = {'','_supra','_infra'};
            if endstatus
                save(fullfile(params.PathName,sprintf('ML-contours_interpolated%s.mat',blade{params.blade+1})),'contours','lines','lines_corrected','status')
            else
                save(fullfile(params.PathName,sprintf('ML-contours%s.mat',blade{params.blade+1})),'contours','lines','status')
            end
            msgbox('The contours have been saved','Save')
        end
        update_plot(F,sizM,Z,lines,contours,status,ui,M)
        
    case 'init'             % initialization at the beginning of the program
        endstatus = false;
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
        if ~isempty(params.FileName(strfind(params.FileName,'_down')+5))
            sample_rate = str2num(params.FileName(strfind(params.FileName,'_down')+5));
        else
            answer = inputdlg('Please give the sampling rate of the stack','Sampling Rate',1,{'1'});
            if isempty(answer)
                return
            end
            sample_rate = str2num(answer{1});
        end
        answer = questdlg('Please give information about the part the stack comprises','Pyramidal Blade','Everything','Only suprapyramidal blade','Only infrapyramidal blade','Everything');
        switch answer
            case 'Everything'
                params.blade = 0;
            case 'Only suprapyramidal blade'
                params.blade = 1;
            case 'Only infrapyramidal blade'
                params.blade = 2;
        end
        if ~isempty(sample_rate) && isnumeric(sample_rate)
            params.sample_rate = sample_rate;
        else
            warndlg('Sample Rate is set to 1')
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
            blade = {'','_supra','_infra'};
            if endstatus 
                 save(fullfile(params.PathName,sprintf('ML-contours_interpolated%s.mat',blade{params.blade+1})),'contours','lines','lines_corrected','status')
            else
                save(fullfile(params.PathName,sprintf('ML-contours%s.mat',blade{params.blade+1})),'contours','lines','status')
            end
%         elseif strcmp(answer,'No')
%             clearvars contours lines status
        elseif strcmp(answer,'Cancel') || isempty(answer)
            return
        end
        if exist('F','var') && ishandle(F)      % close figure
            delete(F)
        end
        
    case 'proceed'
        blade = {'','_supra','_infra'};
        answer = questdlg('Proceed with contour processing and tree-MLyzing?');
        if strcmp(answer,'Yes')
            %             contours = interpz_contours(contours,[],'-ex-w');   % inter- and extrapolate contours to all zplanes
            if ~isempty(lines_corrected)
                if exist(fullfile(params.PathName,sprintf('ML-contours_interpolated%s.mat',blade{params.blade+1})),'file')
                    answer = questdlg('Interpolate z-planes again?');
                else
                    answer = 'Yes';
                end
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
                elseif strcmp(answer,'No')
                    
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
                else
                    return
                end
                contours = makeContours(lines_corrected,1:size(contours,1)); % construct the layer contours from the border lines
                statusx = status;
                lines(1:end-1,:) = lines_corrected([1:(size(lines,1)-1)]*params.Zstep,:);
                lines(end,:) = lines_corrected(end,:);
                status(:) = 6;
                blade = {'','_supra','_infra'};
                save(fullfile(params.PathName,sprintf('ML-contours_interpolated%s.mat',blade{params.blade+1})),'contours','lines','lines_corrected','status')
%                 status = statusx;
                answer = questdlg('Trees already transformed into whole DG?');
                if strcmp(answer,'No')
                    TreeConverter(params.PathName)
                elseif strcmp(answer,'Cancel') || strcmp(answer,'')
                    return
                end
                [FileName,PathName] = uigetfile({'*.mtr','*.mtr - Reconstructed Trees (Treestoolbox)'},'Please open whole DG trees file corresponding to the just segmentated DG',params.PathName);
                if isempty(FileName) || ~ischar(FileName)
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
                        plot_tree(DGtrees{t},DGtrees{t}.R);
                    end
                    for i = 1:5
                        patch(contours{1,i}.Vertices(:,1)*DGtrees{t}.x_scale*params.sample_rate , contours{1,i}.Vertices(:,2)*DGtrees{t}.y_scale*params.sample_rate , repmat(DGtrees{t}.z_scale,[numel(contours{1,i}.Vertices(:,1)) 1]),col{i},'FaceAlpha',0.4);
                        patch(contours{end,i}.Vertices(:,1)*DGtrees{t}.x_scale*params.sample_rate , contours{end,i}.Vertices(:,2)*DGtrees{t}.y_scale*params.sample_rate , repmat(size(contours,1)*DGtrees{t}.z_scale,[numel(contours{end,i}.Vertices(:,1)) 1]),col{i},'FaceAlpha',0.4);
                    end
                    set(gca,'YDir','reverse')
                end
                save_tree({DGtrees},fullfile(PathName,sprintf('%s_MLyzed.mtr',FileName(1:end-4))))    %save MLyzed trees
            end
            endstatus = true;
        end
end

%%
    function update_plot(F,sizM,Z,lines,contours,status,mode,M)     % updates figure
        HP = get(F,'UserData');     % get handles to figure children
        S = get(F,'CurrentAxes');
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
        if endstatus
            set(F,'Name',sprintf('%s - Contours - finished',params.FileName))
        else
            set(F,'Name',sprintf('%s - Contours - in process',params.FileName))
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
            if mode~=1              % also plot all other lines if they exist
                if status(Z) == 6
                    mz = 5;
                    for l = 1:numel(lines{Z,1}.Faces)
                        HP(l+mz+3) = plot(S,cellfun(@(x) x.Vertices(l,1),lines(Z,1:mz)),cellfun(@(x) x.Vertices(l,2),lines(Z,1:mz)),'m-');
                    end
                else
                    mz = status(Z);
                end
                for mo = 1:mz
                    if isfield(lines{Z,mo},'Vertices') && ~isempty(lines{Z,mo}.Vertices)
                        HP(mo+3) = plot(S,lines{Z,mo}.Vertices(:,1),lines{Z,mo}.Vertices(:,2),sprintf('%sx-',col{mo}));
                    end
                end

            end
            if status(Z) ~= 6 && isfield(lines{Z,status(Z)},'Vertices') && ~isempty(lines{Z,status(Z)}.Vertices)
                Point = getPoint(sizM,0);
                if mode ~=1
                    HP(3) = plot(S,lines{Z,status(Z)}.Vertices(:,1),lines{Z,status(Z)}.Vertices(:,2),'r.');   % replot points
                end
                HP(2) = plot(S,[lines{Z,status(Z)}.Vertices(end,1);Point(1)],[lines{Z,status(Z)}.Vertices(end,2);Point(2)],sprintf('%s-',col{status(Z)}));    % plot line to mouse
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
                if all(str2num(answerold{1})==str2num(answer{1}))
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

    function [newGCL newfissura] = make_orthogonal_borders(SGCL,GCL,fissura)
        flag = [false false];
        newGCL(1,:) = GCL(1,:);
        newGCL(size(GCL,1),:) = GCL(end,:);
        newfissura(1,:) = fissura(1,:);
        newfissura(size(fissura,1),:) = fissura(end,:);
        lastp = [1 1];
        curv = sqrt(sum(diff(GCL,2,1).^2,2));    % make pdistr for interp dependent on curvature level
        [nix,ind] = max(curv);
        if abs(ind-size(GCL,1)/2) < size(GCL,1)/20
            newGCL(ceil(size(GCL,1)/2),:) = GCL(ind,:);
            flag(1) = true;
            
            vec = [SGCL(ceil(size(SGCL,1)/2),:) ;SGCL(ceil(size(SGCL,1)/2),:)+(GCL(ind,:)-SGCL(ceil(size(SGCL,1)/2),:))*100];
            [xi,yi,ii] = polyxpoly(vec(:,1),vec(:,2),fissura(:,1),fissura(:,2));   %calculate intersection of SGCL line with norm vector
            if size(ii,1) == 1
                newfissura(ceil(size(fissura,1)/2),:) = fissura(ii(:,2),:);
                flag(2) = true;
            end
        end
% this did not work since curvature of fissura and GCL can be at different
% places
%         curv = sqrt(sum(diff(fissura,2,1).^2,2));    % make pdistr for interp dependent on curvature level
%         [nix,ind] = max(curv);
%         if abs(ind-size(GCL,1)/2) < size(GCL,1)/10
%             newfissura(ceil(size(fissura,1)/2),:) = fissura(ind,:);
%             flag(2) = true;
%         end
        for b = 2:size(SGCL,1)-1
%             flag=false;
            if all(flag) && ceil(size(GCL,1)/2) == b
                flag = [false false];
                lastp(1) = find(all(repmat(newGCL(b,:),[size(GCL,1) 1])==GCL,2));     %refind point of max curv
                lastp(2) = find(all(repmat(newfissura(b,:),[size(fissura,1) 1])==fissura,2));     %refind intersection point
            else
                norm = [1, -(SGCL(b+1,1)-SGCL(b-1,1))/(SGCL(b+1,2)-SGCL(b-1,2))];  %calculate norm vector on SCGL line, das stimmt so!!!
                if isnan(norm(2))
                    errordlg('two points at same place')
                    return
                elseif isinf(norm(2))  % two points have been placed exactly at same y value...
                    norm = [0,1];
                end
                vec = [SGCL(b,:)-norm*10000 ;SGCL(b,:); SGCL(b,:)+norm*10000];
                [xi,yi,ii] = polyxpoly(vec(:,1),vec(:,2),SGCL(:,1),SGCL(:,2));   %calculate intersection of SGCL line with norm vector
                if size(ii,1)>1
                    iii = ii(~all(ii == repmat([2 b],size(ii,1),1),2),:);
                    if iii(1) == 1
                        vec = vec(2:3,:);
                    else
                        vec = vec(1:2,:);
                    end
                end
                [xi,yi,ii] = polyxpoly(vec(:,1),vec(:,2),GCL(:,1),GCL(:,2));   %calculate intersection of GCL line with norm vector
                if all(isempty(ii))
                    if iii(1) == 1
                        vec = [SGCL(b,:)-norm*10000 ;SGCL(b,:)];
                    else
                        vec = [SGCL(b,:); SGCL(b,:)+norm*10000];
                    end
                    [xi,yi,ii] = polyxpoly(vec(:,1),vec(:,2),GCL(:,1),GCL(:,2));   %calculate intersection of GCL line with norm vector
                end
                ind = find(ii(:,2) >= lastp(1));
                if numel(ind)>1  % sollte jetzt gelöst sein
                    [nix, ind] = min(abs(ii(:,2)-b)); % find intersection point which is closer to current site
                end
                if isempty(ind)
                    errordlg('Strange error, change some points and retry')
                    return
%                 elseif ii(ind,2)-lastp(1) > 15       % there might be some problem with finding the correct partner.. if this happens, take the old point and go on
%                     newGCL(b,:) = newGCL(b-1,:);
                else
                    newGCL(b,:) = [xi(ind) yi(ind)];    % this intersection is new point
                    lastp(1) = ii(ind,2);
                end
                
                % %    for debug         figure;plot(newGCL(:,1),newGCL(:,2),'x-'),hold all,plot(SGCL(:,1),SGCL(:,2),'x-'),plot(newfissura(:,1),newfissura(:,2),'-xr')
                %                 figure;plot(newGCL(:,1),newGCL(:,2),'x-'),hold all,plot(SGCL(:,1),SGCL(:,2),'x-'),plot(SGCL(b,1),SGCL(b,2),'-xr')
                %             %             dot([1; -(SGCL(b,1)-SGCL(b-1,1))/(SGCL(b,2)-SGCL(b-1,2))] ,(SGCL(b,:)-SGCL(b-1,:)))
                %             norm1 = [1, -(SGCL(b,1)-SGCL(b-1,1))/(SGCL(b,2)-SGCL(b-1,2))];  %calculate norm vector on SCGL line
                %             if isnan(norm1(2))
                %                 flag = true;
                %             else
                %                 vec = [SGCL(b,:)-norm1*10000 ; SGCL(b,:)+norm1*10000];
                %                 [xi,yi,ii1] = polyxpoly(vec(:,1),vec(:,2),GCL(:,1),GCL(:,2));   %calculate intersection of GCL line with norm vector
                %                 [nix, ind1] = min(abs(ii1(:,2)-b)); % find intersection point which is closer to current site
                % %                 [nix , ind] = min(sqrt(sum((repmat(GCL(b-1,:),numel(xi),1)-[xi,yi]).^2,2)),[],1);
                %                 newGCL(b,:) = [xi(ind1) yi(ind1)];    % this intersection is new point
                %             end
                %             if b~= size(SGCL,1)
                %                 norm2 = [1, -(SGCL(b,1)-SGCL(b+1,1))/(SGCL(b,2)-SGCL(b+1,2))];  %do the same of above again for the norm vector with next GCL point
                %                 if isnan(norm2(2))
                %                     if isnan(norm1(2))
                %                         errordlg('To many points with same coord...')
                %                         return
                %                     end
                %                 else
                %                     vec = [SGCL(b,:)-norm2*10000 ; SGCL(b,:)+norm2*10000];
                %                     [xi,yi,ii2] = polyxpoly(vec(:,1),vec(:,2),GCL(:,1),GCL(:,2));
                %                     [nix, ind2] = min(abs(ii2(:,2)-b)); % find intersection point which is closer to current site
                % %                     [nix , ind] = min(sqrt(sum((repmat(GCL(b-1,:),numel(xi),1)-[xi,yi]).^2,2)),[],1);
                %                     if flag
                %                         newGCL(b,:) = [xi(ind2) yi(ind2)];
                %                     else
                %                         if abs(ii2(ind2,2) - ii1(ind1,2)) > 1   % mean could make a bad coordinate so stay at the GCL
                %                             newGCL(b,:) = GCL(round(mean([ii2(ind2,2),ii1(ind1,2)])),:);
                %                         else
                %                             newGCL(b,:) = mean([newGCL(b,:);[xi(ind2) yi(ind2)]],1);
                %                         end
                %                     end
                %                 end
                %             end
                %now use the point of SGCL and GCL to find orthogonal point at
                %fissura
                vec = [SGCL(b,:) ; SGCL(b,:) + (newGCL(b,:) - SGCL(b,:))*10000];
                [xi,yi,ii] = polyxpoly(vec(:,1),vec(:,2),fissura(:,1),fissura(:,2));   %calculate intersection of fissura line with norm vector
                if size(ii,1)>1
                    [nix,ind] = min(abs(ii(:,2) - lastp(2)));  %look for nearest point
                    xi = xi(ind);
                    yi = yi(ind);
                end
                if ~isempty(xi)
                    newfissura(b,:) = [unique(xi),unique(yi)];
                    lastp(2) = ii(2);
                else
                    newfissura(b,:) = fissura(b,:);
                    lastp(2) = b;
                end
            end
        end
        %both attempts to avoid crossing overs raised problems themselves...
        newfissura = orderline(newfissura);
        newGCL = orderline(newGCL);

        %         [newfissura(:,1),newfissura(:,2)] = poly2cw(newfissura(:,1),newfissura(:,2));
%         newfissura = flipud(newfissura);
%         [newGCL(:,1),newGCL(:,2)] = poly2cw(newGCL(:,1),newGCL(:,2));
%         newGCL = flipud(newGCL);
    end
%%
    function outline = orderline(inline)
%         templine = NaN(size(inline,1),size(inline,2),size(inline,1));
%         for i = 0:size(inline,1)-1
%             templine(:,:,i) = circshift(inline,i);
%         end
%         templine = sqrt(sum(diff(templine,1,3).^2,2));
        outline = zeros(size(inline,1),2);
        outline(1,:) = inline(1,:);
        h = ceil(size(inline,1)/2);
        outline(h,:) = inline(h,:);
        inline(h,:) = [];
%         for i = 1:size(inline,1)-1
%             min(templine(i,1,:))
%         end
    ind = 1;
    for il = 2:size(outline,1)
        if il == h
            continue
        else
            inline(ind,:) = [];
        end
        [nix,ind] = min(sqrt(sum((repmat(outline(il-1,:),[size(inline,1),1])-inline).^2,2)));
        outline(il,:) = inline(ind,:);
    end
    end
end