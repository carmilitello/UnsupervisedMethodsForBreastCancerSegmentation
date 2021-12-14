function [outROI] = regionGrowing_init(imm, mask, opt_T)

[m n] = size(imm); % Dimensions of input image
mask = logical(mask);
reg_size = bwarea(mask); % Number of pixels in region

% % Binary mask erosion to initialize output ROI
% str_el = strel('disk', 2);
% J = imerode(mask, str_el); % Output binary image (ROI)

mask_edge = edge(mask);
[row, col] = find(mask_edge == 1);
% Mean value calculation
imm(~mask) = 0;
% figure; imshow(imm);

reg_mean= nonzeromean(mask.*imm);

% Threshold critetion adaptive calculation
reg_maxdist = abs(opt_T - reg_mean);

% Free memory to store neighbours of the (segmented) region
neg_free = 10000;
%neg_free = Isizes(1) * Isizes(2);
neg_pos = 0;
neg_list = zeros(neg_free,3);

% Distance of the region newest pixel to the region mean
pixdist = 0;

% Neighbor locations (footprint)
neigb = [-1 0; 1 0; 0 -1; 0 1; -1 -1; -1 1; 1 -1; 1 1]; %8-neighborhood

% Sampling of the ROI edge pixels (obtained using the S&M approach)

% Massima grandezza della regione sottoposta a RG (immagine intera)
maxSize= numel(imm);

for b=1:length(row)
    x = row(b); y = col(b);
    
    % Start region-growing until distance between region and posible new pixels
    % become higher than a certain treshold
    while( (pixdist <= reg_maxdist) && (reg_size < maxSize) )
        
        % Add new neighbors pixels
        for j=1:8
            % Calculate the neighbour coordinate (footprint)
            xn = x + neigb(j,1); yn = y + neigb(j,2);
            % Check if neighbour is inside or outside the image
            isInside = ((xn>= 1) && (yn>= 1) && (xn<= m) && (yn<= n));
            
            % Add neighbor if inside and not already part of the segmented area
            if(isInside && (mask(xn,yn)==0))
                neg_pos = neg_pos+1;
                neg_list(neg_pos,:) = [xn yn imm(xn,yn)];
                mask(xn,yn) = 1;
            end
        end
        
        % Add a new block of free memory
        if(neg_pos+10 > neg_free)
            neg_free = neg_free+10000;
            neg_list((neg_pos+1):neg_free,:) = 0;
        end
        
        % Add pixel with intensity nearest to the mean of the region, to the region
        dist = abs(neg_list(1:neg_pos,3) - reg_mean);
        [pixdist, index] = min(dist);
        mask(x,y) = 2;
        reg_size = reg_size+1;
        
        % Calculate the new mean of the region
        reg_mean = (reg_mean*reg_size + neg_list(index,3))/(reg_size+1);
        reg_maxdist = abs(opt_T - reg_mean);
        
        % Save the x and y coordinates of the pixel (for the neighbour add proccess)
        x = neg_list(index,1); y = neg_list(index,2);
        
        % Remove the pixel from the neighbour (check) list
        neg_list(index,:) = neg_list(neg_pos,:);
        neg_pos = neg_pos-1;
        
    end
end

% Return the segmented area as a logical matrix
outROI = mask;

end