function contours = interpz_contours(contours,max_z,options)
if nargin < 1
    [FileName,PathName] = uigetfile('.mat','Choose the contour file to load','ML-contours.mat');
    if FileName ~=0
        load(fullfile(PathName,FileName))
    else
        return
    end
end
if nargin < 2 || isempty(max_z)
    max_z = size(contours,1);
end
if nargin < 3
    options = [];
end

max_c = size(contours,2);

if size(contours,1) < max_z
    contours(max_z,1:max_c)= {[]};
end

if strfind(options,'-w')
    w = waitbar(0,'Interpolating z-planes');
end
for mc=1:max_c
    pairs =[];
    z_cont = zeros(2,1);
    empty = cellfun(@isempty,contours(:,mc));

    ind = find(~empty,1,'first');
    if isempty(ind)
        return
    end
    curr_z = ind -1;
    
    if strfind(options,'-ex')
        contours(1:curr_z,mc) = {[]};
    else
        contours(1:curr_z,mc) = contours(ind,mc);
        fprintf(sprintf('%d-%d\n',1,curr_z));
    end
    
    

    while curr_z <= max_z
        curr_z = curr_z +1;
        fprintf(sprintf('%d\n',curr_z));
        
        if empty(curr_z)
            if isempty(pairs)
                pairs = zeros(max(numel(poly_first.Faces),numel(poly_second.Faces)),2);
                for po = 1:size(poly_first.Vertices,1)
                    dist = sqrt(sum((repmat(poly_first.Vertices(po,:),[size(poly_second.Vertices,1),1])-poly_second.Vertices).^2,2));
                    [nix mindist] = min(dist);
                    pairs(po,:) = [po,mindist];
                end
                remaining_coords = setdiff(1:size(poly_second.Vertices,1),unique(pairs(:,2)));
                for o = 1:numel(remaining_coords)
                    dist = sqrt(sum((repmat(poly_second.Vertices(remaining_coords(o),:),[size(poly_first.Vertices,1),1])-poly_first.Vertices).^2,2));
                    [nix mindist] = min(dist);
                    pairs(po+o,:) = [mindist,remaining_coords(o)];
                end
                pairs = sortrows(pairs);
                fp = sort(pairs(pairs(:,2) == 1 | pairs(:,2) == poly_second.Faces(end),1));
                fp = fp(diff(fp)==0);
                if ~isempty(fp)
                    prob_pairs = pairs(:,1) == fp;
                    partners = pairs(prob_pairs,2);
                    [nix, ind] = max(diff(partners));
                    pairs(prob_pairs,2)=cat(1,partners(ind+1:end),partners(1:ind));
                end
            end
            
            
            contours(curr_z,mc) = interpolate(pairs,poly_first,poly_second,curr_z);
            
        else
            if curr_z == z_cont(2) && isempty(contours{z_cont(1)-1,mc}) && strfind(options,'-ex')
                contours(1:z_cont(1)-1,mc) = interpolate(pairs,poly_first,poly_second,1:z_cont(1)-1,'pchip');
                fprintf(sprintf('%d-%d\n',1,z_cont(1)-1));
                if strfind(options,'-w')
                    waitbar(sum(~cellfun(@isempty,contours(:)))/numel(contours),w);
                end
            end

            ind = find(~empty(curr_z+1:end),1,'first');
            if isempty(ind)
                if sum(~empty) == 1             % then only on contour exists, no interpolation possible
                    contours(:,mc) = contours(curr_z,mc);       % use the only contour for all z_planes
                elseif strfind(options,'-ex')
                    contours(curr_z+1:end,mc) = interpolate(pairs,poly_first,poly_second,curr_z+1:size(contours,1),'pchip');
                else
                    contours(curr_z+1:max_z,mc) = {poly_first};
                end
                fprintf(sprintf('%d-%d\n',curr_z+1,size(contours,1)));
                if strfind(options,'-w')
                    waitbar(sum(~cellfun(@isempty,contours(:)))/numel(contours),w);
                end
                break
            else
                z_cont(1) = curr_z;
                poly_first = contours{curr_z,mc};
                z_cont(2) = ind + curr_z;
                poly_second = contours{z_cont(2),mc};
            end
            pairs =[];
        end
        if strfind(options,'-w')
            waitbar(sum(~cellfun(@isempty,contours(:)))/numel(contours),w);
        end
    end
end
if strfind(options,'-w')
    close(w)
end

if strfind(options,'-s')
    col = {'w','b','g','y','r','k'};
    figure;hold on
    for i=1:size(contours,1)
        for mc = 1:max_c
        this_patch = contours{i,mc};
        this_patch.FaceAlpha = 0.2;
        this_patch.FaceColor = col{mc};
        this_patch.Vertices(:,3)=i;
        patch(this_patch)
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
            this_polygons{z} = this_polygon;
        end
        
    end
end