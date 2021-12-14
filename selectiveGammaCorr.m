function imm = selectiveGammaCorr(imm, alfa)
% Performs a gamma correction only on pixels with an intensity greater than
% the threshold obtained with the Otsu method and sets the others to zero. 
% This is because the lesion is generally lighter and is located to the right
% of the Otsu threshold, while the pizels on the left will not be part of the mask.

    [m n]= size(imm);
    th= graythresh(imm);

    imm= imm(:);
    for i=1:length(imm)
        if imm(i)>= th
            imm(i)= (imm(i))^alfa;
        else
            imm(i)=0;
        end
    end

    imm= reshape(imm, m, n);
    
end
