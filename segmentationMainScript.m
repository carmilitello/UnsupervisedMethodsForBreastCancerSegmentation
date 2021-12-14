
close all;
clear variables;
clc;

disp('Start segmentation');

% list of segmentation algorithm
listAlgorithms = {'sFCM','FCM','kMEANS','SMRG'};

% number of segmentation algorithm
numOfAlgorithms = length(listAlgorithms);
        
% Sets the directory where there are folders with all patients
currFold = uigetdir();

% Get the number of patients (assuming a record for each patient)
allfiles = dir(currFold);

dirFlags = find([allfiles(:).isdir]);
subFolders = allfiles(dirFlags);

% cycle on all patients
for k = 1:length(dirFlags)
    
    % Analizza tutte le fette utili di tutti i pazienti
    if (strcmp(subFolders(k).name, '.') == false && strcmp(subFolders(k).name, '..') == false)
        
        %-----------------------------------------------------------------
        % opening MR images
        paziente = allfiles(k).name;
        
        disp(strcat('  Segmentation of patient:', paziente));

        filenameImages = 'DatasetNIFTI.nii';
        imgpath = fullfile(currFold, paziente, filenameImages);
        
        %read the nifti header
        info = niftiinfo(imgpath);
        
        % Opening the nifti images
        images = niftiread(imgpath);
        
        % image size
        [numRows, numColumns, numSlices] = size(images);

        % In the central slice, the lesion should be more easily visible than the others
        centralSlice = round(numSlices/2);

        % The first and last  slices of each patient do not contain 
        % the lesion and will not be analyzed
        offset = 10;
                
        % select central slice
        selectedImage = images(: ,: , centralSlice);
        selectedImage = imrotate(selectedImage, 90);
        selectedImage = flipdim(selectedImage, 2);

        %------------------------------------------------------------------
        % opening manual masks (it's my ground truth)

        filenameMasksManual = 'DatasetNIFTI_manual_mask.nii';
        manualMaskspath = fullfile(currFold, paziente, filenameMasksManual);

        masksManual = niftiread(manualMaskspath);

        selectedMaskManual = masksManual(: ,: , centralSlice);
        selectedMaskManual = imrotate(selectedMaskManual, 90);
        selectedMaskManual = flipdim(selectedMaskManual, 2); 

        %------------------------------------------------------------------
        % automatic masks opening
        % I need only to follow a variable of the same type.
        % At each iteration (segmentation) then I will overwrite the mask obtained by my algorithm
         
        filenameMaskAuto = 'DatasetNIFTI_mask.nii';
        autoMaskpath = fullfile(currFold, paziente, filenameMaskAuto);

        maskAuto = niftiread(autoMaskpath);

        selectedMaskAuto = maskAuto(: ,: , centralSlice);
        selectedMaskAuto = imrotate(selectedMaskAuto, 90);
        selectedMaskAuto = flipdim(selectedMaskAuto,2); 

        % assign the variable of the same type
        myMasks = maskAuto;
        
        %------------------------------------------------------------------
        % Crop of the area of interest

        im = mat2gray(selectedImage);
        [sizeX, sizeY]= size(im);
        [imcropped, rect] = imcrop(im);
        [m, n] = size(imcropped);
        
        %--------------------------------------------------------------------------
        % cycle on segmentation algorithms
        for algo = 1:numOfAlgorithms

            % Selection of the segmentation algorithm to use
            segmAlgorithm = listAlgorithms{algo};
            
            disp(strcat('    segmentaton with algorithm: ', segmAlgorithm));
            
            %% Segmentation + Postprocessing central slice

            [coord, mask, flag] = segm_centralSlice(imcropped, segmAlgorithm);
            
            % If the lesion is in contact with the crop, double the size
            % of the crop and result in segm_centralSlice

            counter = 0;
            while flag == true

                if 5>= counter
                    rect = [rect(1) rect(2) 2*rect(3) 2*rect(4)];
                    [imcropped, rect] = imcrop(im, rect);
                    [m, n] = size(imcropped);
                    [coord, mask, flag] = segm_centralSlice(imcropped, segmAlgorithm);
                    counter = counter+1;
                else
                    flag = false;
                end
            end

            x0 = coord(1);
            y0 = coord(2); 

            % Return the mask to its correct position
            x_min = round(rect(1))-1;
            y_min = round(rect(2))-1;

            % Display mask (segmentation mask obtained from fcm)
            global_mask = zeros(sizeX, sizeY);
            global_mask(1:m, 1:n) = mask;
            global_mask = imtranslate(global_mask, [x_min, y_min]);

            global_mask_R = imrotate(global_mask, 90);
            global_mask_R = flipdim(global_mask_R, 2);
            
            myMasks(:, :, centralSlice) = global_mask_R;


            %% Segmentation of the slices only with lesion (beyond the central slice)
            
            for selectedSlice = (1+offset):(numSlices-offset)

                if (selectedSlice ~= centralSlice)

                    % selection of the slice under examination from the patient file
                    selectedImage = images(: ,: , selectedSlice);
                    selectedImage = imrotate(selectedImage, 90);
                    selectedImage = flipdim(selectedImage,2);

                    selectedMaskManual = masksManual(: ,: , selectedSlice);
                    selectedMaskManual = imrotate(selectedMaskManual, 90);
                    selectedMaskManual = flipdim(selectedMaskManual,2); 

                    selectedMaskAuto = maskAuto(: ,: , selectedSlice);
                    selectedMaskAuto = imrotate(selectedMaskAuto, 90);
                    selectedMaskAuto = flipdim(selectedMaskAuto,2); 

                    im = mat2gray(selectedImage);
                    [sizeX, sizeY] = size(im);

                    % Crops the same region of interest selected for the center slice
                    [imcropped, rect] = imcrop(im, rect);
                    [m, n] = size(imcropped);

                    %% Segmentation + Postprocessing

                    [mask, flagOther] = segm_otherSlices(imcropped, x0, y0, segmAlgorithm);

                    otherCounter = 0;

                    while flagOther == true
                        if 5>= otherCounter
                            rect = [rect(1) rect(2) 2*rect(3) 2*rect(4)];
                            [imcropped, rect] = imcrop(im, rect);
                            [m, n] = size(imcropped);
                            [mask, flagOther] = segm_otherSlices(imcropped, x0, y0, segmAlgorithm);
                            otherCounter = otherCounter + 1;
                        else
                            flagOther = false;
                        end
                    end

                    %riporta la mask in posizione corretta
                    x_min = round(rect(1))-1;
                    y_min = round(rect(2))-1;

                    %mask visualization
                    global_mask = zeros(sizeX, sizeY);
                    global_mask(1:m, 1:n) = mask;
                    global_mask = imtranslate(global_mask, [x_min, y_min]);

                    global_mask_R = imrotate(global_mask, 90);
                    global_mask_R = flipdim(global_mask_R, 2);
                    
                    myMasks(:, :, selectedSlice) = global_mask_R;

                end
               
            end

            %% Saving results as NIFTI at the end of the cycle on all slices

            filenameMyResult = strcat('DatasetNIFTI_myAlgorithm_', segmAlgorithm, '.nii');
            pathMyResult = fullfile(currFold, paziente, filenameMyResult);
            
            niftiwrite(myMasks, pathMyResult, info);
        
            close all;
        
        end
    end      
end

disp('Finished segmentation');

      