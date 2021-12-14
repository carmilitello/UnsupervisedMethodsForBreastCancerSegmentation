function [coord, mask, flag] = segm_centralSlice(imcropped, segmAlgorithm)
% The function that operates on the central slice proceeds as follows:
% - clustering with FCM, K-means or spatial FCM;
% - delete the components connected with area> areath;
% - eliminates the connected components with eccentricity> .3 or eccentricity> .95;
% - among those left, finds the connected component that maximizes the
% ratio ((MeanIntensity / Eccentricity) * Solidity) and returns the
% centroid.

    % It is initially assumed that the lesion does not touch the edge of the crop
    flag = false;
    
    % Coordinates of the centroid (still zero)
    coord = [0 0];
    [m, n] = size(imcropped);
      
    %% Pre-processing
    
    % Smoothing
    imm = im2double(imgaussfilt(imcropped, .6));
    
    % Centering
    meanValue = nonzeromean(imm);
    imm = imsubtract(imm, meanValue);
     
    % Selective Gamma correction
    imm = selectiveGammaCorr(imm, .6);
    imm = imfill(imm,'holes');

    % Organize data in a row vector (for fcm and Kmeans)
    data = imm(:);
    
    if (strcmp(segmAlgorithm,'FCM'))
        %% FCM
        
        % For the choice of the number of clusters it is evaluated which
        %classification has the lowest partition entropy
        [numClust, flagClust]= optNumClust(data, imm, segmAlgorithm);
        
        options = [2.0 100 1e-5 true];
        [centers, U, obj_fcn] = fcm(data, numClust, options);
        
        % Choice of the highest value of membership functions
        [C, clusters_ID] = max(U);
        
        % mask extraction
        mask = get_mask(centers, clusters_ID, m, n, flagClust);
        
        % image masking
        maskedIm = imcropped .* mask;
         
    elseif (strcmp(segmAlgorithm,'sFCM'))
    %% spatial FCM (sFCM)
    
        [numClust, flagClust]= optNumClust(data, imm, segmAlgorithm);
        
        [U, centers, obj_fcn] = SFCM2D(imm, numClust);
        
        % Choice of the highest value of membership functions
        [C, clusters_ID] = max(U);
        
        % mask extraction
        mask = get_mask(centers, clusters_ID, m, n, flagClust);

        % image masking
        maskedIm = imcropped .* mask;

        
    elseif (strcmp(segmAlgorithm,'kMEANS'))
        %% K-means

        [numClust flagClust]= optNumClust(data, imm, segmAlgorithm);

        [clusters_ID, centers] = kmeans(data, numClust);
        
        % mask extraction
        mask = get_mask(centers, clusters_ID, m, n, flagClust);

        % image masking
        maskedIm = imcropped .* mask;
       
         
    elseif (strcmp(segmAlgorithm,'SMRG'))
        %% Split and merge + region growing
        
        % Since there is no clustering, flagClust = false
        % flagClust= false;        

        % The predicate in entry to splitmerge is based on the comparison
        % between the Otsu threshold of the whole image (opt_T) and the
        % threshold that modifiedOtsu finds for the quadregion in question (m).
        % Splitmerge divides the original image as long as m> opt_T:
        % the exclusion of the darker contributions in the calculation of m
        % allows us to identify a region of greater intensity that probably 
        % represents the lesion.        
        
        [opt_T, EM] = graythresh(imm);
        mask = splitmerge(imm, 2, @new_predicate, opt_T);
        mask = (mask> 0);
        mask = bwareaopen(mask, 5); 
        
        % Threshold change for region growing: given the exclusion of darker
        % contributions, the threshold change allows the growing of the mask
        % previously obtained using a lower threshold value (classic Otsu).
        opt_T = modifiedOtsu(imm);
        
        % Region growing 
        connComp= bwconncomp(mask);     
        numObj= connComp.NumObjects;
        
        for i = 1:numObj
            mask_tmp = zeros(m, n);
            mask_tmp(connComp.PixelIdxList{i}) = 1; 
            J = regionGrowing_init(imm, mask_tmp, opt_T);
            mask = mask + J;
        end

        % If there are more than one connected components, opening with a
        % larger SE eliminates the components that are too small and
        % isolates the weakly connected components
        if numObj>1
            mask= imopen(mask, strel('disk', 2));
        else
            mask = imopen(mask, strel('disk', 1));
        end

        % image masking
        maskedIm = imcropped .* mask;


    end
    
    %% Post- processing

    % Connected components
    cc = bwconncomp(mask);
    extrema = regionprops('table',cc, 'Extrema');
    stats = regionprops('table', cc, imcropped,'Area', 'Centroid', 'Eccentricity', 'MeanIntensity', 'Solidity');

    %% Criterion: calculation of the average area of the connected components
    meanArea= mean(stats.Area, 'all');
    areath= round(3*(meanArea/4));
    underAreas= find(areath>= stats.Area);
    
    % Eliminates connected components smaller than areath pixels
    for i=1: length(underAreas)
        cc.PixelIdxList{1, underAreas(i)}= 0;
        stats(underAreas(i), :)= array2table([0 0 10 0 0]);
    end

    % Intermediate result visualization
    mask = zeros(m, n);
    
    for i=1: length(cc.PixelIdxList)
        if cc.PixelIdxList{i}~= 0
            mask(cc.PixelIdxList{1, i}) = 1;
        end
    end
    
    %% Criterion: eccentricity
    % It is only used if, at this point, there are more than 2 connected 
    % components left for further skimming. Borderline cases have been found
    % in which the lesion had a very elongated shape and its Eccentricity was 
    % very high: since it was the only component left, the next step would have
    % resulted in having 0 connected components to process. 
    % Making the stage below optional avoids this problem.

    if (nnz(stats.Area)>= 2)
        for i=1: length(stats.Eccentricity)
            if ( (.3> stats.Eccentricity(i)) || (stats.Eccentricity(i)> .9) )
                stats(i, :)= array2table([0 0 10 0 0]);
            end
        end

        % Intermediate result visualization
        mask = zeros(m, n);
        candidates = find(stats.Eccentricity~= 10);
        for i=1:length(candidates)
            mask(cc.PixelIdxList{1, candidates(i)}) = 1;
        end
   end
            
    %% Maximize ratio ((MeanIntensity / Eccentricity) * Solidity)
    % Solidity = Area/ConvexArea) makes the lesion search more efficient than Area.
    % Clustering can give rise to very large connected regions other than the lesion:
    % for this reason, the multiplicative factor Area would give an incorrect result.
    % Since the lesions tend to have a convex shape and regular edges,
    % they will have a greater Solidity than other large areas, but with
    % strongly irregular edges: in the latter the Area / ConvexArea ratio
    % will certainly be lower and this region will have a lower ratio ratio.

    ratio = ( (stats.MeanIntensity./stats.Eccentricity).*(stats.Solidity) );
    [maxInt, indmaxInt]= max(ratio, [], 'linear');

    % Intermediate result visualization
    mask = zeros(m, n);
    mask(cc.PixelIdxList{1, indmaxInt})= 1;

    % image masking
    lesion = mask .* imcropped;

    %% ROI edge control
    % If the ROI is in contact with the edge of the crop, it may be that yours
    % shape has been altered: flag == true if there is contact.
    % The flag will be evaluated in main.

    coordinates= cell2mat(extrema.Extrema(indmaxInt));
    minimum = min(coordinates, [],'all');
    maxX= max(coordinates(:, 1), [], 'linear');
    maxY= max(coordinates(:, 2), [], 'linear');

    if ( (1> minimum) || (maxX> n) || (maxY> m))
        flag= true;
    end

    % Stores the centroid of the connected region
    coord(1)= stats.Centroid(indmaxInt, 1);
    coord(2)= stats.Centroid(indmaxInt, 2);

     %% Morphological operations

    mask = imfill(mask, 'holes');
    
end