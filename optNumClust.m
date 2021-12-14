function [numClust, flag]= optNumClust(data, imm, segmAlgorithm)
% optNumClust, based on the chosen segmentation algorithm, evaluates the
% most appropriate number of clusters (in this, between numClust = 2 or 
% numClust = 3): the criterion adopted for the choice is the minimum 
% partition entropy. The definition of vpe was taken from "An automatic 
% fuzzy c-means algorithm for image segmentation", by Y. Li, Y. Shen.

    flag = false;

    % Partition Entropy
    vpe = zeros(1, 2);

    N = length(data);

    if ( (strcmp(segmAlgorithm,'FCM')) || (strcmp(segmAlgorithm,'kMEANS')) )
    %% FCM and kMEANS
    % Since there is no membership function for kMEANS, the number of
    % clusters is estimated as if the chosen algorithm were FCM
    
        % Calculation of the vpe with numClust = 2
        numClust=2;
        [centers, U, obj_fcn] = fcm(data, numClust);
        for i=1: numClust
            for j=1: N
                vpe(1)= vpe(1)+ U(i, j)*log(U(i, j));
            end
            vpe(1)= -vpe(1)/ N;
        end

        % Calculation of the vpe with numClust = 3
        numClust=3;
        [centers , U,  obj_fcn] = fcm(data, numClust);
        for i=1: numClust
            for j=1: N
                vpe(2)= vpe(2)+ U(i, j)*log(U(i, j));
            end
            vpe(2)= -vpe(2)/N;
        end

        [value_vpe numClust_vpe]= min(vpe,[], 'linear');
        numClust= numClust_vpe+1;

    elseif (strcmp(segmAlgorithm,'sFCM'))
    %% sFCM
    
        % Calculation of the vpe with numClust = 2
        numClust=2;
        [U, centers, obj_fcn] = SFCM2D(data, numClust);
        for i=1: numClust
            for j=1: N
                vpe(1)= vpe(1)+ U(i, j)*log(U(i, j));
            end
            vpe(1)= -vpe(1)/ N;
        end

        % Calculation of the vpe with numClust = 3
        numClust=3;
        [U, centers, obj_fcn] = SFCM2D(data, numClust);
        for i=1: numClust
            for j=1: N
                vpe(2)= vpe(2)+ U(i, j)*log(U(i, j));
            end
            vpe(2)= -vpe(2)/N;
        end

        [value_vpe numClust_vpe]= min(vpe,[], 'linear');
        numClust= numClust_vpe+1;

        delta= .01;
        if delta >=abs(vpe(1)- vpe(2))
            numClust=2;
            flag= false;
        end

        if numClust==3
            flag= true;
        end
    
    end
    
end