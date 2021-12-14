function mu = nonzeromean(imm)

    imm = im2double(imm);
    [count, levels] = imhist(imm);

    count = count(2: length(count));
    levels = levels(2: length(levels));
    mu = sum(count.*levels)/sum(count);

end