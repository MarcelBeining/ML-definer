function lines = interpz_lines(lines,max_z,options)
if nargin < 1
    [FileName,PathName] = uigetfile('.mat','Choose the line file to load','ML-lines.mat');
    if FileName ~=0
        load(fullfile(PathName,FileName))
    else
        return
    end
end
if nargin < 2 || isempty(max_z)
    max_z = size(lines,1);
end
if nargin < 3
    options = [];
end

max_c = size(lines,2);

if size(lines,1) < max_z
    lines(max_z,1:max_c)= {[]};
end

if strfind(options,'-w')
    w = waitbar(0,'Interpolating z-planes');
end
for mc=1:max_c
    pairs =[];
    z_cont = zeros(2,1);
    empty = cellfun(@isempty,lines(:,mc));

    ind = find(~empty,1,'first');
    if isempty(ind)
        return
    end
    curr_z = ind -1;
    
    if strfind(options,'-ex')
        lines(1:curr_z,mc) = {[]};
    else
        lines(1:curr_z,mc) = lines(ind,mc);
        if strfind(options,'-c')
            fprintf(sprintf('%d-%d\n',1,curr_z));
        end
    end
    
    

    while curr_z <= max_z
        curr_z = curr_z +1;
        if strfind(options,'-c')
            fprintf(sprintf('%d\n',curr_z));
        end
        if empty(curr_z)
            if isempty(pairs)
                pairs = repmat(poly_first.Faces',[1 2]);%zeros(max(numel(poly_first.Faces),numel(poly_second.Faces)),2);
%                 for po = 1:size(poly_first.Vertices,1)
%                     dist = sqrt(sum((repmat(poly_first.Vertices(po,:),[size(poly_second.Vertices,1),1])-poly_second.Vertices).^2,2));
%                     [nix mindist] = min(dist);
%                     pairs(po,:) = [po,mindist];
%                 end
%                 remaining_coords = setdiff(1:size(poly_second.Vertices,1),unique(pairs(:,2)));
%                 for o = 1:numel(remaining_coords)
%                     dist = sqrt(sum((repmat(poly_second.Vertices(remaining_coords(o),:),[size(poly_first.Vertices,1),1])-poly_first.Vertices).^2,2));
%                     [nix mindist] = min(dist);
%                     pairs(po+o,:) = [mindist,remaining_coords(o)];
%                 end
%                 pairs = sortrows(pairs);
%                 fp = sort(pairs(pairs(:,2) == 1 | pairs(:,2) == poly_second.Faces(end),1));
%                 fp = fp(diff(fp)==0);
%                 if ~isempty(fp)
%                     prob_pairs = pairs(:,1) == fp;
%                     partners = pairs(prob_pairs,2);
%                     [nix, ind] = max(diff(partners));
%                     pairs(prob_pairs,2)=cat(1,partners(ind+1:end),partners(1:ind));
%                 end
            end
            
            
            lines(curr_z,mc) = interpolate(pairs,poly_first,poly_second,curr_z);
            
        else
            if curr_z == z_cont(2) && isempty(lines{z_cont(1)-1,mc}) && strfind(options,'-ex')
                lines(1:z_cont(1)-1,mc) = interpolate(pairs,poly_first,poly_second,1:z_cont(1)-1,'pchip');
                if strfind(options,'-c')
                    fprintf(sprintf('%d-%d\n',1,z_cont(1)-1));
                end
                if strfind(options,'-w')
                    waitbar(sum(~cellfun(@isempty,lines(:)))/numel(lines),w);
                end
            end

            ind = find(~empty(curr_z+1:end),1,'first');
            if isempty(ind)
                if sum(~empty) == 1             % ~then only on contour exists, no interpolation possible
                    lines(:,mc) = lines(curr_z,mc);       % use the only contour for all z_planes
                elseif strfind(options,'-ex')
                    lines(curr_z+1:end,mc) = interpolate(pairs,poly_first,poly_second,curr_z+1:size(lines,1),'pchip');
                else
                    lines(curr_z+1:max_z,mc) = {poly_second};
                end
                if strfind(options,'-c')
                    fprintf(sprintf('%d-%d\n',curr_z+1,size(lines,1)));
                end
                if strfind(options,'-w')
                    waitbar(sum(~cellfun(@isempty,lines(:)))/numel(lines),w);
                end
                break
            else
                z_cont(1) = curr_z;
                poly_first = lines{curr_z,mc};
                z_cont(2) = ind + curr_z;
                poly_second = lines{z_cont(2),mc};
            end
            pairs =[];
        end
        if strfind(options,'-w')
            waitbar(sum(~cellfun(@isempty,lines(:)))/numel(lines),w);
        end
    end
end
if strfind(options,'-w')
    close(w)
end

if 1%strfind(options,'-s')
    col = {'k','y','b','g','r','w'};
    figure;hold on
    for i=1:size(lines,1)
        for mc = 1:max_c
            this_line = lines{i,mc};
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

    function this_polygons = interpolate(pairs,poly_first,poly_second,zs,method)
        if nargin < 5
            method = 'linear';
        end
        this_polygons = cell(numel(zs),1);
        for z = 1:numel(zs)
            curr_z2 = zs(z);
            this_polygon.Vertices = zeros(size(pairs,1),2);
            this_polygon.Faces = 1:size(pairs,1);
            for p = 1:size(pairs,1)
                x = [poly_first.Vertices(pairs(p,1),1),poly_second.Vertices(pairs(p,2),1)];
                y = [poly_first.Vertices(pairs(p,1),2),poly_second.Vertices(pairs(p,2),2)];
                this_polygon.Vertices(p,:) = [interp1([z_cont(1),z_cont(2)],x,curr_z2,method),interp1([z_cont(1),z_cont(2)],y,curr_z2,method)];
            end
% % % % % % %             sum(diff(this_polygon.Vertices(1:end-1,:),1,1).*diff(this_polygon.Vertices(2:end,:),1,1),2)
            this_polygons{z} = this_polygon;
        end
        
    end
end