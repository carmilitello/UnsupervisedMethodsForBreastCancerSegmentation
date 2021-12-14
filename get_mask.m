function mask = get_mask(centers, clusters_ID, m, n, flagClust)
% Create a mask of size mxn:
% - the coordinates of the centroids (centers)
% - the result of the classification (clusters_ID)
% - flagClust (boolean): if flagClust == true-> 3 clusters,
% - if flagClust == false -> 2 clusters

    % Extracts index relating to the cluster of max intensity
    [ROI_center, ROI_index] = max(centers);
    ROI = (clusters_ID == ROI_index);
    mask = reshape(ROI, m, n);
    mask = imopen(mask, strel('disk', 1));

    if flagClust
        % Search for the second cluster with higher intensity
        centers(ROI_index) = 0;
        [ROI_center2, ROI_index2] = max(centers);
        ROI2 = (clusters_ID == ROI_index2);
        mask2 = reshape(ROI2, m, n);
        % OR logico tra la mask del primo e secondo cluster
        % (per includere anche regione necrotica)
        mask = or(mask, mask2);
        % Morphological opening per eliminare le parti non rilevanti della
        % mask (poorly connected regions)
        mask= imopen(mask, strel('disk', 2)); 
    end
        
end