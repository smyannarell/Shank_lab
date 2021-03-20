clear

%setup variables
path_errors = 0;
files_array = ["testcode"];
error_paths = {};
array_location = 1;
x = 0;
y = 0;
output_array = {};
output_row = 2;
output_collumn = 1;
next_location = 1;
per_reporter_matrix = struct();
per_reporter_matrix.name = ["test"];
errors_list = struct();
errors_list.location = ["no errors found"];
error_list.type = ["no errors found"];
number_of_errors = 0;
location_in_files_array = 1;
iii = 1;
figure_location = 1;
output_path = "D:\Outputs\";

%error codes
error1 = "path invalid";
error2 = "unknown microscopy channel";
error3 = "same data twice";
error4 = "some strains are missing";
error5 = "resolution error";

%this is used for troubleshooting further down the pipe. If it's ever >0
%check naming convention.

%import strain number to gene reporter list 
[~,~,Strain_information] = xlsread('C:\Users\Eric\Documents\Work\Shank Lab\Sarah Y\2019- Shank lab strains with pairwise reporters for flow cytometry.xlsx');

%Select folder that contains all the images to be segmented
image_directory = dir(uigetdir('.','Where are your images located?'));

%Index folder
for specific_image=3:length(image_directory)
    if contains(image_directory(specific_image).name, "DIC")
        x=1;
    elseif contains(image_directory(specific_image).name, "YFP")
        x=2;
    elseif contains(image_directory(specific_image).name, "mTurq")
        x=3;
    end
    
    if contains(image_directory(specific_image).name, "8X")
        x=x+3;
    end
    
    if contains(image_directory(specific_image).name, "SY")
        if contains(image_directory(specific_image).name, '_')
            SY_number = extractBetween(image_directory(specific_image).name,'S','_','Boundaries','exclusive');
            date_done = extractBefore(image_directory(specific_image).name,'S');
        elseif contains(image_directory(specific_image).name, '-');
            SY_number = extractBetween(image_directory(specific_image).name,'S','-','Boundaries','exclusive');
            date_done = extractBefore(image_directory(specific_image).name,'S');
        else
            number_of_errors = number_of_errors + 1;
            errors_list.location(number_of_errors,1) = sprintf(image_directory(specific_image).name);
        end
            
        if any(any(contains(files_array,SY_number)&contains(files_array,date_done)))
           [current_SY_location,trashvariable] = find(contains(files_array,SY_number)&contains(files_array,date_done), 1,'first');
           files_array{current_SY_location,x} = sprintf(image_directory(specific_image).name); 
        else
            files_array{next_location,x} = sprintf(image_directory(specific_image).name);
            next_location = next_location+1;
            stupidmatlabstep = cellfun('isempty',files_array);
            files_array(stupidmatlabstep) = cellstr("Is Missing");
        end
    else
        number_of_errors = number_of_errors + 1;
        errors_list.location(number_of_errors,1) = sprintf(image_directory(specific_image).name);
        errors_list.type(number_of_errors,1) = sprintf(error1);
    end
end


%%
%Identify Gene/reporter combo
for location_in_files_array=1:length(files_array)
%for location_in_files_array = 10
    if any(any(contains(files_array(location_in_files_array,:),"Is Missing")))
        output_row = output_row+1;
        number_of_errors = number_of_errors + 1;
        errors_list.location(number_of_errors,1) = sprintf("Row number" + output_row + " is missing an output");
        errors_list.type(number_of_errors,1) = sprintf(error4);
        
    else
        if contains(files_array(location_in_files_array,1),'_')
            Current_SY_number = extractBetween(files_array(location_in_files_array,1),'S','_','Boundaries','exclusive');
        elseif contains(files_array(location_in_files_array,1),'-')
            Current_SY_number = extractBetween(files_array(location_in_files_array,1),'S','-','Boundaries','exclusive');
        end
        
        Location_of_SY_number = contains(Strain_information(:,2),Current_SY_number);
        YFP_name = extractBetween(Strain_information(Location_of_SY_number,3),'E::','-Ypet','Boundaries','exclusive');
        mTurq_name = extractBetween(Strain_information(Location_of_SY_number,3),'A::','-mTurq','Boundaries','exclusive');

        %Import Images
        DIC = rgb2gray(imread(files_array(location_in_files_array,1)));
        YFP = rgb2gray(imread(files_array(location_in_files_array,2))); 
        mTurq = rgb2gray(imread(files_array(location_in_files_array,3)));

        %report if images are not the same resolution
        if size(DIC)~=size(YFP)|size(DIC)~=size(mTurq)|size(YFP)~=size(mTurq)
            number_of_errors = number_of_errors + 1;
            errors_list.location(number_of_errors,1) = sprintf(Current_SY_number);
            errors_list.type(number_of_errors,1) = sprintf(error5);
        end

        %Auto-Threshold/Segement Images
        %DIC image
        gaussDIC = imgaussfilt(DIC, 1);
        grayDIC = mat2gray(gaussDIC);
        DIC_threshold = graythresh(grayDIC);
        segDIC = grayDIC>DIC_threshold;
        
        %YFP image
        gaussYFP = imgaussfilt(YFP, 1);
        grayYFP = mat2gray(gaussYFP);
        YFP_threshold = graythresh(grayYFP);
        segYFP = grayYFP>YFP_threshold;
        
        %mTurq image
        gaussmTurq = imgaussfilt(mTurq, 1);
        graymTurq = mat2gray(gaussmTurq);
        mTurq_threshold = graythresh(graymTurq);
        segmTurq = graymTurq>mTurq_threshold;
        
        %Fill in holes in DIC image to find center
        filled_DIC = imfill(segDIC, 'holes');


        %Identify total pixels in colony
        DICarea = sum(sum(segDIC));

        %Identify total pixels in colony that express reporter
        YFParea = sum(sum(and(segYFP,segDIC)));
        mTurqarea = sum(sum(and(segmTurq,segDIC)));
        coexpressedarea = sum(sum(and(and(segYFP,segmTurq),segDIC)));

        %find center of DIC
        %[r,c] = find(segDIC == 1);
        %max_y = max(r);
        %min_y = min(r);
        %max_x = max(c);
        %min_x = min(c);
        %center = [round((((max_x - min_x)/2)+min_x),0), round((((max_y - min_y)/2)+min_y),0)];
        image_statistics_Area = regionprops(filled_DIC, 'Area');
        image_statistics = regionprops(filled_DIC);
        all_image_centers = struct2array(image_statistics_Area);
        [~,main_colony_center] = max(all_image_centers);
        center = [round(image_statistics(main_colony_center).Centroid)];
        
        %calculate euclidean distance between each point in the image
        [cartsian_x,cartsian_y] = meshgrid(1:size(DIC,2),1:size(DIC,1));
        eucledian_distance = sqrt((center(1)-cartsian_x).^2+(center(2)-cartsian_y).^2);
        eucledian_distance = round(eucledian_distance,0);

        %Look for overlap between Reporter A and DCI
        percentYFP = YFParea/DICarea*100;

        %Look for overlap between Reporter B and DCI
        percentmTurq = mTurqarea/DICarea*100;

        %Look for coexpression
        percentcoexpressed = coexpressedarea/DICarea*100;

        %Quantify average expression in each region
        averageYFP = mean(YFP(and(segYFP,segDIC)));
        averagemTurq = mean(mTurq(and(segmTurq,segDIC)));
        averagemTurqincoexpressed = mean(mTurq(segmTurq&segDIC&segYFP));
        averageYFPincoexpressed = mean(YFP(segmTurq&segDIC&segYFP));

        %Average reading over distance from center of colony
        linearized_euclidian_distance=eucledian_distance(:);
        linearized_mTurq=mTurq(:);
        linearized_YFP=YFP(:);
        [Unique_Euclidian_distance,~,row_location_of_linearized_euclidian_distance]=unique(linearized_euclidian_distance);
        mean_signal_over_distance_mTurq=accumarray(row_location_of_linearized_euclidian_distance,linearized_mTurq,[],@mean);
        mean_signal_over_distance_YFP=accumarray(row_location_of_linearized_euclidian_distance,linearized_YFP,[],@mean);

        vi = 1;
        standard_deviation_list_mTurq = [];
        standard_deviation_list_YFP = [];
        Confidence_interval_mTurq = [];
        Confidence_interval_YFP = [];
        for vi = 1:length(Unique_Euclidian_distance)
            Rows_of_a_given_distance = row_location_of_linearized_euclidian_distance == vi;
            number_of_pixes_at_given_length = sum(Rows_of_a_given_distance);
            target_T_score = tinv([0.025 0.975],number_of_pixes_at_given_length - 1);

            standard_deviation_in_mTurq = std(double(linearized_mTurq(Rows_of_a_given_distance)));
            standard_deviation_list_mTurq(vi,:) = [vi standard_deviation_in_mTurq];
            Standard_error_mTurq = standard_deviation_list_mTurq(vi,2)/sqrt(number_of_pixes_at_given_length);
            Confidence_interval_mTurq(vi,:) = Standard_error_mTurq * target_T_score;

            standard_deviation_in_YFP = std(double(linearized_YFP(Rows_of_a_given_distance)));
            standard_deviation_list_YFP(vi,:) = [vi standard_deviation_in_YFP];
            Standard_error_YFP = standard_deviation_list_mTurq(vi,2)/sqrt(number_of_pixes_at_given_length);
            Confidence_interval_YFP(vi,:) = Standard_error_YFP * target_T_score;
        end

        %Image tracing
        boundaries_list = bwboundaries(filled_DIC);
        all_boundaries = boundaries_list{:};
        eucledian_distance_of_boundaries = sqrt((center(1)-all_boundaries(:,1)).^2+(center(2)-all_boundaries(:,2)).^2);
        eucledian_distance_of_boundaries = round(eucledian_distance_of_boundaries,0);
        average_pixels_to_edge_of_colony = mean(eucledian_distance_of_boundaries);
        
        %calculate percent to edge of colony
        percent_to_edge_of_colony = (Unique_Euclidian_distance/average_pixels_to_edge_of_colony)*100;
        
        %graph average signal over distance from center of colony
        if any(any(contains(files_array(location_in_files_array,:),"Is Missing")))

        else
            h(figure_location) = figure('Name',"S" + string(Current_SY_number));
            tiledlayout(2,1)
            nexttile
            plot(percent_to_edge_of_colony,mean_signal_over_distance_mTurq);
            hold on
            
            %Below lines will plot standard deviation
            %plot(Unique_Euclidian_distance,mean_signal_over_distance_mTurq+standard_deviation_list_mTurq(:,2));
            %plot(Unique_Euclidian_distance,mean_signal_over_distance_mTurq-standard_deviation_list_mTurq(:,2));
            
            %Below lines will calculate 95% confidence interval
            plot(percent_to_edge_of_colony,mean_signal_over_distance_mTurq + Confidence_interval_mTurq);
            
            title(mTurq_name + " -- mTurq");
            xlabel("percent distance to edge of colony")
            ylabel("average intensity of pixels at distance")
            hold off
            nexttile
            plot(percent_to_edge_of_colony,mean_signal_over_distance_YFP);
            hold on
            
            %Below lines will plot standard deviation
            %plot(Unique_Euclidian_distance,mean_signal_over_distance_YFP+standard_deviation_list_YFP(:,2));
            %plot(Unique_Euclidian_distance,mean_signal_over_distance_YFP-standard_deviation_list_YFP(:,2));
            
            %below lines will calculate 95% confidence interval
            plot(percent_to_edge_of_colony,mean_signal_over_distance_YFP + Confidence_interval_YFP);
            
            title(YFP_name + " -- Ypet");
            xlabel("percent distance to edge of colony");
            ylabel("average intensity of pixels at distance");
            hold off
            sgtitle("S" + string(Current_SY_number))
            set(gcf, 'Position', [100, 100, 600, 600])
            figure_location = figure_location+1;
        end


        end
end

%Save all figures for later use
savefig(h, output_path+'AllDistanceFigures_24hr.fig')
close(h)

fname = "output_24hr.csv";
%fname = fullfile(PathName,filename);

writecell(output_array, output_path + fname);