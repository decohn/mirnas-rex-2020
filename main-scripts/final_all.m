%import normal data
a = fileDatastore('**/*.txt', 'ReadFcn', @importdata)
filenames = a.Files
numberfiles = length(filenames)
normaloutput = array2table(zeros(1881,71))
normalnameoutput = strings(1,71)
%%
for k = 1: numberfiles
    fprintf('Now reading file %s\n', filenames{k}); 
    normaldata = readtable(filenames{k});
    output = normaldata(:,3); %take only the reads per million
    normaloutput(:,k) = output 
    [filepath, name] = fileparts(filenames{k})
    normalnameoutput(:,k) = name
end
%%
load('normaloutput')
writetable(normaloutput, 'normaloutput.csv') %don't need to do this yet
load('normalnameoutput')
writematrix(normalnameoutput, 'normalnameoutput.csv')
%%
%added .txt to the end and then transposed manually
input1 = readtable('normalnameoutput.csv')
normalnames = table2array(input1)
input2 = readtable('normal_sample_sheet.txt')
normalsamplesheet = table2array(input2)
%%
%compare coloumns and replace 
[x,y] = ismember(normalnames(:,1),normalsamplesheet(:,2));
out = normalnames;
out(x,2) = normalsamplesheet(y(x),6)
%%
%manually checked in excel using conditionally formating for duplicates =
%none
finalnormal = readtable('finalnormaloutput.csv', 'ReadVariableNames', true)
%%
%importing tumor data: split into 5 smaller sections as taking too long to
%run in one big file
a = fileDatastore('**/*.txt', 'ReadFcn', @importdata)
filenames = a.Files
numberfiles = length(filenames)
tumoroutput1 = array2table(zeros(1881,108))
tumornameoutput1 = strings(108,1)
%%
for k = 1:numberfiles
    fprintf('Now reading file %s\n', filenames{k}); 
    tumordata1 = readtable(filenames{k});
    output = tumordata1(:,3); %take only the reads per million
    tumoroutput1(:,k) = output 
    [filepath, name] = fileparts(filenames{k})
    tumornameoutput1(:,k) = name
end
%%
a = fileDatastore('**/*.txt', 'ReadFcn', @importdata)
filenames = a.Files
numberfiles = length(filenames)
tumoroutput2 = array2table(zeros(1881,119))
tumornameoutput2 = strings(1,119)
%%
for k = 1:numberfiles
    fprintf('Now reading file %s\n', filenames{k}); 
    tumordata2 = readtable(filenames{k});
    output = tumordata2(:,3); %take only the reads per million
    tumoroutput2(:,k) = output 
    [filepath, name] = fileparts(filenames{k})
    tumornameoutput2(:,k) = name
end
%%
a = fileDatastore('**/*.txt', 'ReadFcn', @importdata)
filenames = a.Files
numberfiles = length(filenames)
tumoroutput3 = array2table(zeros(1881,129))
tumornameoutput3 = strings(1,129)
%%
for k = 1:numberfiles
    fprintf('Now reading file %s\n', filenames{k}); 
    tumordata3 = readtable(filenames{k});
    output = tumordata3(:,3); %take only the reads per million
    tumoroutput3(:,k) = output 
    [filepath, name] = fileparts(filenames{k})
    tumornameoutput3(:,k) = name
end
%%
a = fileDatastore('**/*.txt', 'ReadFcn', @importdata)
filenames = a.Files
numberfiles = length(filenames)
tumoroutput4 = array2table(zeros(1881,105))
tumornameoutput4 = strings(1,105)
%%
for k = 1: numberfiles
    fprintf('Now reading file %s\n', filenames{k}); 
    tumordata4 = readtable(filenames{k});
    output = tumordata4(:,3); %take only the reads per million
    tumoroutput4(:,k) = output 
    [filepath, name] = fileparts(filenames{k})
    tumornameoutput4(:,k) = name
end
%%
a = fileDatastore('**/*.txt', 'ReadFcn', @importdata)
filenames = a.Files
numberfiles = length(filenames)
tumoroutput5 = array2table(zeros(1881,83))
tumornameoutput5 = strings(1,83)
%%
for k = 1: numberfiles
    fprintf('Now reading file %s\n', filenames{k}); 
    tumordata5 = readtable(filenames{k});
    output = tumordata5(:,3); %take only the reads per million
    tumoroutput5(:,k) = output 
    [filepath, name] = fileparts(filenames{k})
    tumornameoutput5(:,k) = name
end
%%
%save name files to add .txt and transpose in excel
writematrix(tumornameoutput1, 'tno1.xlsx')
writematrix(tumornameoutput2, 'tno2.xlsx')
writematrix(tumornameoutput3, 'tno3.xlsx')
writematrix(tumornameoutput4, 'tno4.xlsx')
writematrix(tumornameoutput5, 'tno5.xlsx')
%%
%load all files
input1 = readtable('tno1.xlsx', 'ReadVariableNames', false)
name1 = table2array(input1)
input2 = readtable('tno2.xlsx', 'ReadVariableNames', false)
name2 = table2array(input2)
input3 = readtable('tno3.xlsx', 'ReadVariableNames', false)
name3 = table2array(input3)
input4 = readtable('tno4.xlsx', 'ReadVariableNames', false)
name4 = table2array(input4)
input5 = readtable('tno5.xlsx', 'ReadVariableNames', false)
name5 = table2array(input5)
input6 = readtable('gdc_sample_sheet.2020-12-20.txt')
samples = table2array(input6)
%%
%finding patient IDs
[x,y] = ismember(name1(:,1),samples(:,2));
out1 = name1;
out1(x,2) = samples(y(x),6);
[x,y] = ismember(name2(:,1),samples(:,2));
out2 = name2;
out2(x,2) = samples(y(x),6);
[x,y] = ismember(name3(:,1),samples(:,2));
out3 = name3;
out3(x,2) = samples(y(x),6);
[x,y] = ismember(name4(:,1),samples(:,2));
out4 = name4;
out4(x,2) = samples(y(x),6);
[x,y] = ismember(name5(:,1),samples(:,2));
out5 = name5;
out5(x,2) = samples(y(x),6)
%%
load('to1')
writetable(tumoroutput1, 'to1.csv')
load('to2')
writetable(tumoroutput2, 'to2.csv')
load('to3')
writetable(tumoroutput3, 'to3.csv')
load('to4')
writetable(tumoroutput4, 'to4.csv')
load('to5')
writetable(tumoroutput5, 'to5.csv')
%%
finaltumor = readtable('finaltumoroutput.csv', 'ReadVariableNames', true)
%%
%load all csv output
normal = readtable('finalnormaloutput.csv', 'ReadVariableNames', true)
tumor = readtable('finaltumoroutput.csv', 'ReadVariableNames', true)
%%
%selecting only paired data
logical = ismember(tumor.Properties.VariableNames, normal.Properties.VariableNames)
pairedtumor = tumor(:, logical == 1)
%%
%sort the data tables alphabetically
nsort = sort(normal.Properties.VariableNames(2:end));
sortednormal = [normal(:,1) normal(:,nsort)]
tsort = sort(pairedtumor.Properties.VariableNames(2:end));
sortedtumor = [pairedtumor(:,1) pairedtumor(:,tsort)]
%%
%discard miRNAs expressed average <1 in both normal and tumor samples
mnormal = mean(sortednormal{:,2:end},2)
mtumor = mean(sortedtumor{:,2:end},2)
lcombined = mnormal <1 & mtumor <1
%%
%only keep where the above logical conditions are false
finalnormal = sortednormal(lcombined == 0,:)
finaltumor = sortedtumor(lcombined == 0,:)
load('finalnormal')
writetable(finalnormal, 'finalnormal.xlsx')
load('finaltumor')
writetable(finaltumor, 'finaltumor.xlsx')
%%
%new script and workspace starting here
input1 = readtable('finaltumor.xlsx', 'ReadVariableNames', true)
input2 = readtable('finalnormal.xlsx', 'ReadVariableNames', true)
tumor = table2array(input1(:,2:end))
normal = table2array(input2(:,2:end))
Nrows = size(tumor,1)
h = zeros(Nrows,1)
p_value = zeros(Nrows,1)
for k = 1:Nrows
  [h(k), p_value(k)] = ttest(tumor(k,:),normal(k,:))
end
%%
%fold change
mtumor = mean(input1{:,2:end},2)
mnormal = mean(input2{:,2:end},2)
N_rows = size(mtumor,1)
Fold_Change = zeros(N_rows,1)
for o = 1:N_rows
  Fold_Change(o) = mnormal(o,:)/mtumor(o,:)
end
%%
%combining everything into a table
x = [Fold_Change p_value]
y = array2table(x)
names = input1(:,1)
z = [names y]
result = sortrows(z,3)
%% 
%now apply Benjamini-Hochberg correction
rank = [1:1:387]'
rank1 = array2table(rank)
ttest1 = [result rank1]
%%
FDR1 = (0.05*ttest1.rank)/387
FDR2 = array2table(FDR1)
ttestfinal = [result FDR2]
%%
%manual inspection of ttestfinal table indicates that critical pvalue is
%2.86e-27
check = ttestfinal.x2 <= 2.861609440892827e-27
final = ttestfinal(check == 1,:)
%%
%extract only columns needed and export as .csv
needed = final(:,[1:3])
needed.Properties.VariableNames = {'Name' 'Fold_Change' 'p_value'}
writetable(needed, 'export.csv')