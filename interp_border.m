% FUNCTION interp_border interpolates the contour "line" either by
% distributing a defined amount of points in equidistance (mode 1),in a defined
% distance (mode 2) or dependent on a point distribution histogram. Additionally, smoothing can be applied to the contour.
% if no pdistr is given, a gaussian distribution is shown

function outline = interp_border(line,x,mode,options,filtsize,pdistr)

if nargin < 4 || isempty(options)
    options = '';
end
if mode == 3 && (nargin < 6 || isempty(pdistr))
    pdistr = -100:100;
    pdistr = 10/(20*sqrt(2*pi)) * exp(-0.5*(pdistr/(2*50)).^2);
end

outline = [];
if size(line,1) < 2  %too low number of points
    return
end


d_line = diff(line,1,1);    % calc vectors from one to next point
lengths = sqrt(sum(d_line.^2,2));   % calc lengths of these vectors
d_line = d_line./repmat(lengths,[1 size(d_line,2)]);    % norm the vectors
glength = sum(lengths); % calc length of whole contour
if mode == 1
    if x < 2    %too low number of points
        return 
    end
    pnumber = x;
    ilength = glength/(pnumber-1);  % calc distance between interpolated points
elseif mode == 2
    if x > glength  % too large intervall
        return
    end
    ilength = x;
    pnumber = round(glength/ilength)+1; % calc number of intervall points 
    % CAUTION: End point of interpolated line may change drastically if
    % intervall length is huge!
elseif mode == 3
    pnumber = x;
    if numel(pdistr) ~= pnumber     %interpolate distribution if it does not fit to number of points
        rat = (numel(pdistr)-1)/(pnumber-1);
        pdistr = interp1(1:numel(pdistr),pdistr,1:rat:numel(pdistr),'pchip');
    end
    pdistr = -pdistr + max(pdistr); % invert the distribution so that huge values will give narrow points
    pdistr = pdistr / sum(pdistr);  % norm it again
else
    return
end
if mode == 3
    cum_ilength = cumsum(pdistr*glength);
else
    cum_ilength = cumsum(repmat(ilength,[pnumber-1,1]));    %calc cumulative sum of these distances
end

cum_ilength(end) = glength;     % avoid calc-error when cumulating ilengths

cum_lengths = cumsum([0;lengths]);      %calc distances of original points to first point (along contour)

outline = zeros(pnumber,size(line,2));  % predefine output array
outline(1,:) = line(1,:);   %first point is original point

for p = 1:pnumber-1
    ind = find(cum_ilength(p)>cum_lengths,1,'last');    % find the original point which is closest to the current distance
    outline(p+1,:) = line(ind,:) + d_line(ind,:)*(cum_ilength(p)-cum_lengths(ind));   %calc position of interpolated point by using the closest original point and the vector to next point
end


if strfind(options,'smooth')
    %     filtsize = 0.05;
    if round(pnumber*filtsize) ~= pnumber*filtsize
        kern = rem(floor(pnumber*filtsize),2) * floor(pnumber*filtsize) + rem(ceil(pnumber*filtsize),2) * ceil(pnumber*filtsize); % rounds to the odd number
    else
        kern = pnumber*filtsize + ~rem(pnumber*filtsize ,2) ;
    end
    tmp = outline([1,end],:);
    outline(:,1) = convn (padarray(outline(:,1),(kern-1)/2,'replicate'), ones (1, kern)' / kern, 'valid');
    outline(:,2) = convn (padarray(outline(:,2),(kern-1)/2,'replicate'), ones (1, kern)' / kern, 'valid');
    outline([1,end],:) = tmp;   % accounts for smoothing error at edges
end