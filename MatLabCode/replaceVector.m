function [Vnew] = replaceVector(Vorig,Vrep,indrep)
% Function that takes vector Vorig and replaces the elements at indices indrep
% with the vectors contained in the cell Vrep

% Initialize new vector
untouched_ind_old = setdiff([1:length(Vorig)]',indrep);
sizeVnew = length(untouched_ind_old);
for j=1:size(indrep,1)
    sizeVnew = sizeVnew + size(Vrep{j},1);
end
Vnew = zeros(sizeVnew,1);
% Fill in Vnew
indOLD = 1;
indNEW = 1;
for i=1:size(indrep,1)
    % Fill in untouched elements up to index of replacement
    length_untouched = indrep(i)-indOLD;
    Vnew(indNEW:indNEW+length_untouched-1) = Vorig(indOLD:indOLD+length_untouched-1);
    indNEW = indNEW + length_untouched;
    indOLD = indOLD + length_untouched;
    % Fill in new elements
    length_touched = length(Vrep{i});
    Vnew(indNEW:indNEW+length_touched-1) = Vrep{i};
    indNEW = indNEW + length_touched;
    indOLD = indOLD + 1;
end
% Fill in last untouched chunk (if applicable)
if indOLD <= length(Vorig)
    Vnew(indNEW:end) = Vorig(indOLD:end);
end

end