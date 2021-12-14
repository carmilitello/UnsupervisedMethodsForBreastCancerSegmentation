function [MAD, MAXD] = distance_metrics(reference, automated)
%DISTANCE_METRICS Function that calculates and returns the average distance
%(average segmentation error) and the maximum distance (maximum segmentation
% error) between the boundaries (edges) between an image segmented by the
% automatic method and a reference image (manual segmentation performed by 
% a rediologist).

    b_ref = bwboundaries(reference);
    b_ref = b_ref{1};
    N = length(b_ref);

    b_auto = bwboundaries(automated);
    b_auto = b_auto{1};
    K = length(b_auto);

    % Distance vector
    dists = zeros(K,1);


    for i = 1:K
        d_min = 300;
        for j = 1:N
            d = sqrt( (b_auto(i,1) - b_ref(j,1) )^2 + (b_auto(i,2) - b_ref(j,2) )^2 );

            if (d < d_min)
                d_min = d;
            end

        end
        dists(i) = d_min;
    end

    % Mean Absolute Distance
    MAD = mean(dists);

    % Max Absolute Distance
    MAXD = max(dists);

end

