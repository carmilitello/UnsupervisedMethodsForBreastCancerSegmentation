function flag = new_predicate(imm, th)

    % Predicate based on the mean of the quadregion in question
    m = modifiedOtsu(imm);

    % Binarization according to the Otsu method: the part to the right of the 
    %threshold will certainly contain the gray levels affecting the lesion
    flag = (m > th) ;
end
