contours =[];
contours{5}.Faces = [1 2 3 4 5];
contours{5}.Vertices = [2 4;10 12;20 24;6 20;4 10];
contours{10}.Faces = [1 2 3 4 5];
contours{10}.Vertices = [3 10;2 5;1 2;5 6;10 12];
max_z = 20;


%%
contours{max_contours}= [];
pairs =[];

empty = cellfun(@isempty,contours);

ind = find(~empty,1,'first');
if isempty(ind)
    return
end
curr_z = ind -1;
contours(1:curr_z) = contours(ind);

while curr_z <= max_z
    curr_z = curr_z +1;
    fprintf(sprintf('%d\n',curr_z));
    
    if empty(curr_z)
        if isempty(pairs)
            pairs = zeros(max(numel(poly_first.Faces),numel(poly_second.Faces)),2);
            for p = 1:size(poly_first.Vertices,1)
                dist = sqrt(sum((repmat(poly_first.Vertices(p,:),[size(poly_second.Vertices,1),1])-poly_second.Vertices).^2,2));
                [nix mindist] = min(dist);
                pairs(p,:) = [p,mindist];
            end
            remaining_coords = setdiff(1:size(poly_second.Vertices,1),unique(pairs(:,2)));
            for o = 1:numel(remaining_coords)
                dist = sqrt(sum((repmat(poly_second.Vertices(remaining_coords(o),:),[size(poly_first.Vertices,1),1])-poly_first.Vertices).^2,2));
                [nix mindist] = min(dist);
                pairs(p+o,:) = [mindist,remaining_coords(o)];
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
            this_polygon.Faces = 1:size(pairs,1);
        end
        
        
        this_polygon.Vertices = zeros(size(pairs,1),2);
        for p = 1:size(pairs,1)
            x = [poly_first.Vertices(pairs(p,1),1),poly_second.Vertices(pairs(p,2),1)];
            y = [poly_first.Vertices(pairs(p,1),2),poly_second.Vertices(pairs(p,2),2)];
            this_polygon.Vertices(p,:) = [interp1([z_cont(1),z_cont(2)],x,curr_z),interp1([z_cont(1),z_cont(2)],y,curr_z)];
        end
        contours{curr_z} = this_polygon;
        
    else
        poly_first = contours{curr_z};
        z_cont(1) = curr_z;
        ind = find(~empty(curr_z+1:end),1,'first');
        if isempty(ind)
            contours(curr_z+1:max_z) = {poly_first};
            break
        else
            z_cont(2) = ind + curr_z;
            poly_second = contours{z_cont(2)};
        end
        pairs =[];
    end
end

%%
figure;hold on
for i=1:numel(contours)
    this_patch = contours{i};
    this_patch.FaceAlpha = 0.3;
    this_patch.FaceColor = [1 0 1];
    this_patch.Vertices(:,3)=i;
   patch(this_patch)
end
   