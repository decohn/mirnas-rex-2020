%% Setting up the data



GDC_sample_sheet_tumour = readtable("../wetransfer-a1573d/tumour_sample_file.txt");

collate_files_into_single_folder("../wetransfer-a1573d/Tumours");

sample_rpm = readtable("../wetransfer-a1573d/Tumours/00a99edf-57d6-49f2-b810-d9de193c2881.mirbase21.isoforms.quantification.txt");

list_Sample_ID_tumour = GDC_sample_sheet_tumour.CaseID;



length(GDC_sample_sheet_tumour.SampleID)
%544
%number of samples (columns)
length(sample_rpm.reads_per_million_miRNA_mapped)
%2089
%number of mirnas (rows) 



list_tumour_files = dir("../wetransfer-a1573d/Tumours");
list_tumour_files(3) = [];
list_tumour_files(2) = [];
list_tumour_files(1) = [];


%% importing the sample sheet, averaging, removing duplicates 

rpm_matrix = zeros(2089, 544);

for current_file = 1:length(list_tumour_files)
    %something weird in file, skip first 3 lines
    current_file_data = readtable("../wetransfer-a1573d/Tumours/" +list_tumour_files(current_file).name);
        for gdc_sample_sheet_row = 1:length(GDC_sample_sheet_tumour.FileName)
            if strcmp(list_tumour_files(current_file).name, GDC_sample_sheet_tumour.FileName(gdc_sample_sheet_row)) == 1
                sample_id = list_Sample_ID_tumour(gdc_sample_sheet_row);
                
                
                rpm_matrix(:, gdc_sample_sheet_row) = current_file_data.reads_per_million_miRNA_mapped;  % accesses the element in the first row and second column
                 
             
            end 
         end 
end


for x = 1:544
    checked_ids = zeros(615); %setting to empty matrix
                                %once for loop is run, we will alter
                                %checked_ids, so duplicates have 1
    duplicate_id = list_Sample_ID_tumour(x); %accessing the list 
    denominator = 0; %setting count to 0 
    reads_per_million_sum = zeros(2089); % setting another empty matrix (for 
                                        %for sums of reads per million for
                                        %duplicate columns
                                    
    
    %making a sum of reads per millions for duplicated data
    
    for y = 1:544
        if checked_ids(y) == 0 && strcmp(duplicate_id, list_Sample_ID_tumour(y))
            reads_per_million_sum(:,1) = reads_per_million_sum(:,1) + rpm_matrix(:,y);
            % X(:,y) is the matching column sharing the same sample id but
            %w/ a different tumour sample
            denominator = denominator + 1;
            % count's the numbers of duplicates 
        end
    end
    
    reads_per_million_sum = reads_per_million_sum / denominator;
    
    %creating new column for averaged 
    
    for y = 1:544
        if checked_ids(y) == 0 && strcmp(duplicate_id, list_Sample_ID_tumour(y))
            rpm_matrix(:,y) = reads_per_million_sum(:,1);
            %if id hasn't been checked, and if duplicate id =
            %list_sample_id, then replace column data (changing value of x)
            %w/ averaged data
            checked_ids(y) = 1;
            %once this is done, mark the column as checked w/ a 1. 
        end
    end
    
end


%% removing duplicates in list
unique_list_Sample_ID = unique(list_Sample_ID_tumour, "stable");

rpm_unique = unique(rpm_matrix.','rows', "stable").';
rpm_table = array2table(rpm_unique, "VariableNames",  string(unique_list_Sample_ID));

rpm_table.mirna_name = rand(2089,1);
rpm_table.mirna_name = sample_rpm.miRNA_ID;

mirna_names = rpm_table.mirna_name;
%% clinical

clinical = readtable("../clinical.cases_selection.2021-01-06/clinical.txt");
clinical = [clinical(:,2) clinical(:,10) clinical(:,16) clinical(:,48)];

%davids method
%%if alive 
clinical.days = zeros(1074,1);
patient_alive = ismember(clinical.vital_status, "Alive"); 


days_for_living_patients = clinical.days_to_last_follow_up(patient_alive);
%days to last follow up for only patients who are still ALIVE

days_for_living_patients = cellstr(num2str(days_for_living_patients));
days_for_living_patients = cellfun(@str2double, days_for_living_patients);

clinical.days(patient_alive) = days_for_living_patients;
%filters based on patient_alive list 


%%%if dead
patient_dead = ismember(clinical.vital_status, "Dead"); 

days_for_dead_patients = clinical.days_to_death(patient_dead);
%days to death only exists for patients who are dead 

days_for_dead_patients = cellfun(@str2double, days_for_dead_patients);

clinical.days(patient_dead) = days_for_dead_patients;

clinical = unique(clinical);

%clinical.Properties.RowNames =  clinical.case_submitter_id


%% making sure only data that exists for both clinical + miRNA are included
string_unique_list_Sample_ID = string(unique_list_Sample_ID);

char_unique_list_Sample_ID = convertStringsToChars(string_unique_list_Sample_ID);

clinical_and_mirna = intersect(char_unique_list_Sample_ID, clinical.case_submitter_id);

rpm_table = rpm_table(:,clinical_and_mirna);

%rpm_table.mirna_name = rand(2089,1);
rpm_table.mirna_name = sample_rpm.miRNA_ID;

for clinical_row = length(clinical.case_submitter_id):-1:1
    if  ismember(clinical.case_submitter_id(clinical_row), clinical_and_mirna) == 0
        clinical(clinical_row,:) = [];
    end
    
end 
        
for clinical_row = length(clinical.case_submitter_id):-1:1
    if clinical_row>1 && strcmp(clinical.case_submitter_id(clinical_row),clinical.case_submitter_id(clinical_row - 1))
       clinical(clinical_row,:) = [];
    end
end 





%% brians mirnas - removing mirnas that aren't from brians list


brians_mirnas = readtable("../mirnas-rex-2020-main 2/output-files/de-mirnas.csv");

brians_mirna_names = brians_mirnas.Name;
filtered_mirna_names = rpm_table.mirna_name;

shared_mirna_names = intersect(filtered_mirna_names, brians_mirna_names);

for mirna_row = length(rpm_table.mirna_name):-1:1
    if  ismember(rpm_table.mirna_name(mirna_row), shared_mirna_names) == 0
        rpm_table(mirna_row,:) = [];
    end
    
end 

%% performing survival analysis
set(groot,'defaultFigureVisible','off')
p_values_from_MatSurv = zeros(210, 1);

for rpm_column = 1:length(rpm_table.mirna_name)
    MatSurvData = MatSurv(clinical.days, clinical.vital_status, table2array(rpm_table(rpm_column, 1:end-1)), 'MedianLess', false);
    if ~isempty(MatSurvData)
        p_values_from_MatSurv(rpm_column) = MatSurvData;
    end
end

saved_data = p_values_from_MatSurv;

%% performing Benjamini-Hochberg correction
saved_data = array2table(saved_data);
saved_data.mirna_name = brians_mirnas.Name;


saved_data.result = fdr_bh(saved_data.saved_data, 0.05, 'pdep', 'yes');

 
 for final_miRNA = length(saved_data.saved_data):-1:1
     if saved_data.result(final_miRNA) == 0
         saved_data(final_miRNA, :) = [];
     end
     
 end 



saved_data = [saved_data(:,1) saved_data(:,2)];

saved_data.Properties.VariableNames(1) = {'p_value'};


