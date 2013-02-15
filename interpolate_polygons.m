contours{1}.Faces = [1 2 3 4 5];
contours{1}.Vertices = [1 2;5 6;10 12;3 10;2 5];
contours{10}.Faces = [1 2 3 4 5];
contours{10}.Vertices = [1 2;5 6;10 12;3 10;2 5];

% contours = cell(100,1);

%%
max_contours = 100;
curr_contours=0;
not_empty = ~cellfun(@isempty,contours);
while curr_contours <= max_contours
    curr_contours = curr_contours +1;
    fprintf(sprintf('%d\n',curr_contours));
%     if curr_contours > numel(contours)
%         break
%     end
    if isempty(contours{curr_contours})
        pairs =[];
        this_polygon= [];
        contour1 = find(not_empty(1:curr_contours),1,'last');
        contour2 = find(not_empty(curr_contours:end),1,'first')+curr_contours-1;
        
        if isempty(contour1) 
            contour1 = find(not_empty,1,'first');
            curr_contours = contour1 -1;
            contours(1:curr_contours) = contours(contour1);
            not_empty(1:curr_contours) = true;
            continue
        elseif isempty(contour2)
            contours(curr_contours:max_contours) = contours(curr_contours-1);
            not_empty(curr_contours:max_contours) = true;
            break
        end
        
        poly_first = contours{contour1};
        poly_second = contours{contour2};
        
        for p = 1:sicontourse(poly_first,1)
            dist = sqrt(sum((repmat(poly_first(p,:),[sicontourse(poly_second,1),1])-poly_second).^2));
            [nix mindist] = min(dist);
            pairs(end+1,:) = [p,mindist];
        end
        remaining_coords = setdiff([1:sicontourse(poly_second,1)],unique(pairs(:,2)));
        for p = 1:numel(remaining_coords)
            dist = sqrt(sum((repmat(poly_second(remaining_coords(p),:),[sicontourse(poly_first,1),1])-poly_first).^2));
            [nix mindist] = min(dist);
            pairs(end+1,:) = [mindist,remaining_coords(p)];
        end
        
        for p = 1:sicontourse(pairs,1)
            x = [poly_first(pairs(p,1),1),poly_second(pairs(p,2),1)];
            y = [poly_first(pairs(p,1),2),poly_second(pairs(p,2),2)];
            this_polygon(end+1,:) = [interp1([contour1,contour2],x,curr_contours),interp1([contour1,contour2],y,curr_contours)];
        end
        contours{curr_contours} = this_polygon;
        not_empty(curr_contours) = true;
    end
end

if curr_contours ~= max_contours
    errordlg('Möp')
end