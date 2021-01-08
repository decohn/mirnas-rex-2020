%import normal data
a = fileDatastore('**/*.txt', 'ReadFcn', @importdata)
filenames = a.Files
numberfiles = length(filenames)
normaloutput = array2table(zeros(2089,71))
normalnameoutput = strings(1,71)
%%
for k = 1: numberfiles
    fprintf('Now reading file %s\n', filenames{k}); 
    normaldata = readtable(filenames{k});
    output = normaldata(:,2); %take only the reads per million
    normaloutput(:,k) = output 
    [filepath, name] = fileparts(filenames{k})
    normalnameoutput(:,k) = name
end
%%
writematrix(normalnameoutput, 'normalnameoutput.xlsx')
%% 
%added.txt to strings and transposed manually
normalnames = readtable('normalnameoutput.xlsx','ReadVariableNames', false) 
%%
%compare coloumns and replace 
normalsamplesheet = readtable('normal_sample_file.txt', 'ReadVariableNames', true)
%%
[x,y] = ismember(normalnames.Var1,normalsamplesheet.FileName);
out = normalnames;
out(x,2) = normalsamplesheet(y(x),6)
%%
%copy and paste names into the normal output
writetable(normaloutput, 'finalnormal.xlsx')
%%
%importing tumor data
a = fileDatastore('**/*.txt', 'ReadFcn', @importdata)
filenames = a.Files
numberfiles = length(filenames)
tumoroutput1 = array2table(zeros(2089,125))
tumornameoutput1 = strings(1,125)
%%
for k = 1:numberfiles
    fprintf('Now reading file %s\n', filenames{k}); 
    tumordata1 = readtable(filenames{k});
    output = tumordata1(:,2); %take only the reads per million
    tumoroutput1(:,k) = output 
    [filepath, name] = fileparts(filenames{k})
    tumornameoutput1(:,k) = name
end
%%
%importing tumor data
a = fileDatastore('**/*.txt', 'ReadFcn', @importdata)
filenames = a.Files
numberfiles = length(filenames)
tumoroutput2 = array2table(zeros(2089,93))
tumornameoutput2 = strings(1,93)
%%
for k = 1:numberfiles
    fprintf('Now reading file %s\n', filenames{k}); 
    tumordata2 = readtable(filenames{k});
    output = tumordata2(:,2); %take only the reads per million
    tumoroutput2(:,k) = output 
    [filepath, name] = fileparts(filenames{k})
    tumornameoutput2(:,k) = name
end
%%
%importing tumor data
a = fileDatastore('**/*.txt', 'ReadFcn', @importdata)
filenames = a.Files
numberfiles = length(filenames)
tumoroutput3 = array2table(zeros(2089,85))
tumornameoutput3 = strings(1,85)
%%
for k = 1:numberfiles
    fprintf('Now reading file %s\n', filenames{k}); 
    tumordata3 = readtable(filenames{k});
    output = tumordata3(:,2); %take only the reads per million
    tumoroutput3(:,k) = output 
    [filepath, name] = fileparts(filenames{k})
    tumornameoutput3(:,k) = name
end
%%
%importing tumor data
a = fileDatastore('**/*.txt', 'ReadFcn', @importdata)
filenames = a.Files
numberfiles = length(filenames)
tumoroutput4 = array2table(zeros(2089,108))
tumornameoutput4 = strings(1,108)
%%
for k = 1:numberfiles
    fprintf('Now reading file %s\n', filenames{k}); 
    tumordata4 = readtable(filenames{k});
    output = tumordata4(:,2); %take only the reads per million
    tumoroutput4(:,k) = output 
    [filepath, name] = fileparts(filenames{k})
    tumornameoutput4(:,k) = name
end
%%
%importing tumor data
a = fileDatastore('**/*.txt', 'ReadFcn', @importdata)
filenames = a.Files
numberfiles = length(filenames)
tumoroutput5 = array2table(zeros(2089,133))
tumornameoutput5 = strings(1,133)
%%
for k = 1:numberfiles
    fprintf('Now reading file %s\n', filenames{k}); 
    tumordata5 = readtable(filenames{k});
    output = tumordata5(:,2); %take only the reads per million
    tumoroutput5(:,k) = output 
    [filepath, name] = fileparts(filenames{k})
    tumornameoutput5(:,k) = name
end
%%
%added .txt and transposed
writematrix(tumornameoutput1, 'tno1.xlsx')
writematrix(tumornameoutput2, 'tno2.xlsx')
writematrix(tumornameoutput3, 'tno3.xlsx')
writematrix(tumornameoutput4, 'tno4.xlsx')
writematrix(tumornameoutput5, 'tno5.xlsx')
%%
tno1 = readtable('tno1', 'ReadVariableNames', false)
tno2 = readtable('tno2', 'ReadVariableNames', false)
tno3 = readtable('tno3', 'ReadVariableNames', false)
tno4 = readtable('tno4', 'ReadVariableNames', false)
tno5 = readtable('tno5', 'ReadVariableNames', false)
%%
tumorsamplesheet = readtable('tumour_sample_file.txt')
%%
[x,y] = ismember(tno1.Var1,tumorsamplesheet.FileName);
out1 = tno1;
out1(x,2) = tumorsamplesheet(y(x),6)
[x,y] = ismember(tno2.Var1,tumorsamplesheet.FileName);
out2 = tno2;
out2(x,2) = tumorsamplesheet(y(x),6)
[x,y] = ismember(tno3.Var1,tumorsamplesheet.FileName);
out3 = tno3;
out3(x,2) = tumorsamplesheet(y(x),6)
[x,y] = ismember(tno4.Var1,tumorsamplesheet.FileName);
out4 = tno4;
out4(x,2) = tumorsamplesheet(y(x),6)
[x,y] = ismember(tno5.Var1,tumorsamplesheet.FileName);
out5 = tno5;
out5(x,2) = tumorsamplesheet(y(x),6)
%%
%manually get rid of duplicates and add patient IDs
writetable(tumoroutput1, 'to1.xlsx')
writetable(tumoroutput2, 'to2.xlsx')
writetable(tumoroutput3, 'to3.xlsx')
writetable(tumoroutput4, 'to4.xlsx')
writetable(tumoroutput5, 'to5.xlsx')
%%
finaltumor = readtable('to1.xlsx', 'ReadVariableNames', true)
finalnormal = readtable('finalnormal.xlsx', 'ReadVariableNames', true)
%%
%selecting only paired data
logical = ismember(finaltumor.Properties.VariableNames, finalnormal.Properties.VariableNames)
pairedtumor = finaltumor(:, logical == 1)
%%
%sort tables alphabetically
nsort = sort(finalnormal.Properties.VariableNames(2:end));
sortednormal = [finalnormal(:,1) finalnormal(:,nsort)]
tsort = sort(pairedtumor.Properties.VariableNames(2:end));
sortedtumor = [pairedtumor(:,1) pairedtumor(:,tsort)]
%%
%discard miRNAs expressed average <1 in both normal and tumor samples
mnormal = mean(sortednormal{:,2:end},2)
mtumor = mean(sortedtumor{:,2:end},2)
lcombined = mnormal <1 & mtumor <1
%%
%only keep where the above logical conditions are false
normal = sortednormal(lcombined == 0,:)
tumor = sortedtumor(lcombined == 0,:)
%%
writetable(normal, 'finalnormal.xlsx')
writetable(tumor, 'finaltumor.xlsx')
%%
%new script and workspace starting here and do ttest
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
%combining everything into a table and sort rows in ascending order
x = [Fold_Change p_value]
y = array2table(x)
names = input1(:,1)
z = [names y]
result = sortrows(z,3)
%%
%now apply Benjamini-Hochberg correction
rank = [1:1:495]'
rank1 = array2table(rank)
ttest1 = [result rank1]
%%
FDR1 = (0.05*ttest1.rank)/495
FDR2 = array2table(FDR1)
ttestfinal = [result FDR2]
%%
%manual inspection of ttestfinal table indicates that critical pvalue is
%in row 404
check = ttestfinal(1:404,:)
%%
%apply fold change filter for final results
fcfilter = check.x1 >= 2.0 | check.x1 <= 0.5
results = check(fcfilter == 1, :)
%%
%extract only coloumns needed and export as .csv
needed = results(:,[1:3])
needed.Properties.VariableNames = {'Name' 'Fold_Change' 'p_value'}
writetable(needed, 'de-mirnas.csv')