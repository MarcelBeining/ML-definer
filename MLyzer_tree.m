function trees = MLyzer_tree(trees,contours,options)

if nargin < 3
    options = [];
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
        coord_plane = find(trees{t}.Z >= (Z-1) * z_scale & trees{t}.Z < Z * z_scale);
        if ~isempty(coord_plane)
            isGCL = inpolygon(trees{t}.X(coord_plane),trees{t}.Y(coord_plane),contours{Z,1}.Vertices(:,1) * x_scale,contours{Z,1}.Vertices(:,2) * y_scale);
            isIML = inpolygon(trees{t}.X(coord_plane),trees{t}.Y(coord_plane),contours{Z,2}.Vertices(:,1) * x_scale,contours{Z,2}.Vertices(:,2) * y_scale);
            isMML = inpolygon(trees{t}.X(coord_plane),trees{t}.Y(coord_plane),contours{Z,3}.Vertices(:,1) * x_scale,contours{Z,3}.Vertices(:,2) * y_scale);
            isOML = inpolygon(trees{t}.X(coord_plane),trees{t}.Y(coord_plane),contours{Z,4}.Vertices(:,1) * x_scale,contours{Z,4}.Vertices(:,2) * y_scale);
            trees{t}.R(coord_plane(isGCL)) = 2;
            trees{t}.R(coord_plane(isIML)) = 3;
            trees{t}.R(coord_plane(isMML)) = 4;
            trees{t}.R(coord_plane(isOML)) = 5;
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
    trees{t}.rnames(2:5) = {'GCL','IML','MML','OML'};
end