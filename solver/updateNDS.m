function [indexlist, add_index] = updateNDS(nds, new_point) 
%updateNDS Updates a nondominated set nds given a new_point

size_nds = size(nds,2);
if isempty(nds)
    % If there is no point in nds add new_point
    add_index = 1;
    indexlist = [true];
elseif any(all(bsxfun(@ge,new_point,nds)))
    % If newPoint is dominated by an entry from ListPNS, it will not be
    % added
    add_index = false;
    indexlist = true(1,size_nds);
else
    % Add newPoint to ListPNS
    add_index = size_nds+1;
    % Delete all entries of ListPNS that are dominated by newPoint
    indexlist = [~all(bsxfun(@le,new_point,nds)),true];
end
end