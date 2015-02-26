function trees = MLyzer_tree(trees,contours,sample_rate,options)

if nargin < 4
    options = [];
end
if nargin < 3 || isempty(sample_rate)
    if isfield(contours{1,1},'sample_rate')
        sample_rate = contours{1,1}.sample_rate;
    else
        sample_rate = 1;
    end
end
structflag = false;
if isstruct(trees)
    trees = {trees};
    structflag = true;
end

x_scale = unique(cellfun(@(x) single(x.x_scale),trees));
y_scale = unique(cellfun(@(x) single(x.y_scale),trees));
z_scale = unique(cellfun(@(x) single(x.z_scale),trees));

if numel(x_scale)>1 || numel(y_scale)>1 || numel(z_scale)>1
    errordlg(sprintf('Some trees have different z_scale than others'))
    return
end

if strfind(options,'-w')
    w = waitbar(0,'MLyzing each tree z-plane');
end

for Z = 1:size(contours,1)
    for t = 1: numel(trees)
        if Z ==1
            coord_plane = find(trees{t}.Z < Z * z_scale);
        elseif Z == size(contours,1)
            coord_plane = find(trees{t}.Z >= (Z-1) * z_scale & trees{t}.Z < Z * z_scale);
        else
            coord_plane = find(trees{t}.Z >= (Z-1) * z_scale);
        end
        if ~isempty(coord_plane)
            isSGCL = inpolygon(trees{t}.X(coord_plane),trees{t}.Y(coord_plane),contours{Z,1}.Vertices(:,1) * sample_rate * x_scale,contours{Z,1}.Vertices(:,2) * sample_rate * y_scale);
            isGCL = inpolygon(trees{t}.X(coord_plane),trees{t}.Y(coord_plane),contours{Z,2}.Vertices(:,1) * sample_rate * x_scale,contours{Z,2}.Vertices(:,2) * sample_rate * y_scale);
            isIML = inpolygon(trees{t}.X(coord_plane),trees{t}.Y(coord_plane),contours{Z,3}.Vertices(:,1) * sample_rate * x_scale,contours{Z,3}.Vertices(:,2) * sample_rate * y_scale);
            isMML = inpolygon(trees{t}.X(coord_plane),trees{t}.Y(coord_plane),contours{Z,4}.Vertices(:,1) * sample_rate * x_scale,contours{Z,4}.Vertices(:,2) * sample_rate * y_scale);
            isOML = inpolygon(trees{t}.X(coord_plane),trees{t}.Y(coord_plane),contours{Z,5}.Vertices(:,1) * sample_rate * x_scale,contours{Z,5}.Vertices(:,2) * sample_rate * y_scale);
            trees{t}.R(coord_plane(isSGCL)) = 1;
            trees{t}.R(coord_plane(isGCL)) = 2;
            trees{t}.R(coord_plane(isIML)) = 3;
            trees{t}.R(coord_plane(isMML)) = 4;
            trees{t}.R(coord_plane(isOML)) = 5;
            trees{t}.R(coord_plane(~isSGCL & ~isGCL & ~isIML & ~isMML & ~isOML)) = 6;
        end
    end
    if strfind(options,'-w')
        waitbar(Z/size(contours,1),w);
    end
end
if strfind(options,'-w')
    close(w)
end

for t = 1: numel(trees)
    trees{t}.rnames(1:6) = {'SGCL','GCL','IML','MML','OML','outside'};
end

if structflag
    trees = trees{1};
end