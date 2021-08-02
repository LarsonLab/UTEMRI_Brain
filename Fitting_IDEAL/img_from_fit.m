function [map] = img_from_fit(result, mask)

map = zeros(size(mask));
map(mask) = result;

end