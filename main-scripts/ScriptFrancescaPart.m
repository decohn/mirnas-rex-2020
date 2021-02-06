clc;
clear;

%%%%%%%%%%%%%%%%%%---------- miRNA ------------%%%%%%%%%%%%%%%%%%%%%%%%

%Table with miRNA cases, case ID and file names,...
miRNAinitialtable = readtable ( "tumour_sample_file.txt");

%Pull all files in openable folder
collate_files_into_single_folder('Tumours');

%Pull all filenames out of initial tablle into vector
miRNAfilenames=miRNAinitialtable(2:end,2);
miRNAtotcases=size(miRNAfilenames); % determine total number of miRNA cases
miRNAtotcases=miRNAtotcases(1);

%Pull out only first file name to get number of miRNA number
miRNAfirstfilename = string(miRNAfilenames{1,1});
miRNAfirstdata = readtable(miRNAfirstfilename); % pull data for first miRNA filename
miRNAID=miRNAfirstdata(:, 1); %Get all miRNA ID's from first file
totmiRNA=size(miRNAID); %Get total number fo miRNA's
totmiRNA=totmiRNA(1);

miRNAcasedatamatrix=zeros(totmiRNA,miRNAtotcases);

 list=ls('Tumours');
% 
 for i=1:miRNAtotcases
    casefilename = string(miRNAfilenames{i,1});
     casedata = readtable(casefilename);
     casemildata=casedata(:,2);
     casemildata=table2array(casemildata);
     miRNAcasedatamatrix(:,i)=casemildata;
 end

% miRNAcasedatamatrix is the matrix with all the data from the everyone part
%%%%%%%%%%%%%%%%%%%--------------------------%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%
%%%Start of my section

%Need to do step 1 to 3 when I have the data and can input in the website

% % A is the fake output that was already provided
% A = readtable('de-mirnas.csv');
% 
% totmiRNA2=size(A);
% totmiRNA2=totmiRNA2(1);
% simplifiedFinalMatrix=zeros(totmiRNA2,miRNAtotcases); 
% L = table2cell(casedata);
% for k=1:totmiRNA2
%     miRNAname=A(k,1);
%     C=table2cell(miRNAname);
%     for kk=1:totmiRNA
% 
%         U= strcmpi(C ,L(kk, 1));
%         if U==1
%             simplifiedFinalMatrix(k,:)= miRNAcasedatamatrix(kk,:); 
%         end
%     end
% end

%%%%%%%%%%%%%%%%%%------------ mRNA -----------------%%%%%%%%%%%%%%%%%%

%dataFrancesca is the sample sheet that I downloaded for my section
mRNAinitialtabledata = readtable("Perfect Francesca Data Sheet");

collate_files_into_single_folder('Francesca TCGA-KIRC Data');
% The new data is GoodData Francesca TCGA-KIRC Data--because Francesca
% TCGA-KIRC Data the files are zipped

mRNAfilenames=mRNAinitialtabledata(2:end,8);
mRNAtotcases=size(mRNAfilenames);
mRNAtotcases=mRNAtotcases(1);
mRNAfirstfilename = string(mRNAfilenames{1,1});
mRNAfirstdata = readtable (mRNAfirstfilename);
mRNAID=mRNAfirstdata(:, 1);
totmRNA=size(mRNAID);
totmRNA=totmRNA(1);


mRNAcasedatamatrix=zeros(totmRNA,mRNAtotcases);

listFrancesca=ls('GoodData Francesca TCGA-KIRC Data'); 

for i=1:mRNAtotcases
    casefilenameFrancesca = string(mRNAfilenames{i,1});
    casedataFrancesca = readtable(casefilenameFrancesca);
    casemildataFrancesca=casedataFrancesca(:,2);
    casemildataFrancesca=table2array(casemildataFrancesca);
    mRNAcasedatamatrix(:,i)=casemildataFrancesca;
end

%mRNAcasedatamatrix is the matrix with all the data for my part (mRNA
%matrix) 
%%%%%%%%%%%%%%%%%%----------------------------%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%----- Selection only differentiated miRNA's ------%%%%%%%%%%%
differentiatedmiRNAsbrian = readtable('DifferentiallyExpressedmiRNAsBrian.txt');

totdifferentiatedmiRNA=size(differentiatedmiRNAsbrian);
totdifferentiatedmiRNA=totdifferentiatedmiRNA(1);

miRNAcasedatamatrixdifferentiated=zeros(totdifferentiatedmiRNA,miRNAtotcases); 
miRNA_ID = table2cell(casedata);
for k2=1:totdifferentiatedmiRNA
    miRNAname2=differentiatedmiRNAsbrian(k2,1);
    C2=table2cell(miRNAname2);
    for kk2=1:totmiRNA

        O= strcmpi(C2 ,miRNA_ID(kk2, 1));
        if O==1
            miRNAcasedatamatrixdifferentiated(k2,:)= miRNAcasedatamatrix(kk2,:); 
        end
    end
end
%miRNAcasedatamatrixdifferentiated is the new miRNA final matrix
%%%%%%%%%%%%%%%--------------------%%%%%%%%%%%%%%%%%%%%%%%%%



%% PART I AM CURRENTLY WORKING ON 

% mirDIP performed and the resulting file with Data is
%This file contains the miRNA and their Gene Symbol target
mRNATargetbymiRNADiff = readtable("miRNATArget.xlsx");

%e!Ensembl was used to determine all Gene Symbols for all Gene Stable ID
% All the files FPKM and the mRNA matrix is formed by Gene Stable ID
EnsemblData= readtable("Ensembl Data .txt"); 


%% %% STEP-BY-STEP

%% 1. Use EnsemblData to convert the Gene Stable IDs in mRNAID into gene names.

% Everything appearing after a dot in a Gene Stable ID is a version number
% that isn't constant between different releases of Ensembl. We need to get
% rid of these, so that we can ensure consistency between the GSIDs in the
% TCGA data and the GSIDs in the file you downloaded from Ensembl recently.
mRNAIDCell = table2cell(mRNAID);
mRNAID2 = regexprep(mRNAIDCell, '\..*', '');

EnsemblData.GeneStableIDVersion = regexprep(EnsemblData.GeneStableID, '\..*', '');
EnsemblData2= table2cell(EnsemblData);

%Ensembl Data cell is the tabloe which has Ensembl Data-Gene Stable ID,
%Gene Name and Gene Stable ID
%mRNAID2 has Gene Stable ID


% Try to find a way to do the conversion. This is a bit tricky; I recommend
% using ismember(), and especially its second output (LocB). If you can't
% figure this out, let me know tomorrow and I'll send you code for it. Try
% starting with 

 [~, locs] = ismember(mRNAID2, EnsemblData.GeneStableIDVersion);
 
 matchingmRNA = EnsemblData.GeneName(locs(locs ~= 0));
 mRNAcasedatamatrix2 = mRNAcasedatamatrix(locs ~= 0, :);


%mRNAcasedatamatrix2 is the matrix with mRNA that is present on Ensembl and
%previous matrix
        




%%
 miRNA_samples_names = miRNAinitialtable.TCGA_A3_3387;
 mRNA_sample_names = mRNAinitialtabledata.Var11;
 shared_sample_names = intersect(miRNA_samples_names, mRNA_sample_names);
 [cc,a]=ismember(miRNA_samples_names, shared_sample_names);
 [sorted_vara, sorting_indicesa] = sort(a);
 sorted_vara=unique(sorted_vara);
%  miRNAinitialtable2= miRNAinitialtable(sorted_vara,:);
 miRNAinitialtable= miRNAinitialtable(ismember(miRNA_samples_names, mRNA_sample_names), :);
 mRNAinitialtabledata = mRNAinitialtabledata(ismember(mRNA_sample_names, shared_sample_names), :);
 
%  
 miRNA_ID = miRNAinitialtable.TCGA_A3_3387;
 mRNA_ID = mRNAinitialtabledata.Var11;
[shared_caseID, ImiRNA, ImRNA] = intersect(miRNA_ID, mRNA_ID,'rows');
miRNAinitialtable=miRNAinitialtable(ImiRNA,:);
mRNAinitialtabledata=mRNAinitialtabledata(ImRNA,:);

miRNAfilenames=miRNAinitialtable(2:end,2);
miRNAtotcases=size(miRNAfilenames); % determine total number of miRNA cases
miRNAtotcases=miRNAtotcases(1);

%Pull out only first file name to get number of miRNA number
miRNAfirstfilename = string(miRNAfilenames{1,1});
miRNAfirstdata = readtable(miRNAfirstfilename); % pull data for first miRNA filename
miRNAID=miRNAfirstdata(:, 1); %Get all miRNA ID's from first file
totmiRNA=size(miRNAID); %Get total number fo miRNA's
totmiRNA=totmiRNA(1);

miRNAcasedatamatrix=zeros(totmiRNA,miRNAtotcases);

list=ls('Tumours');

for i=1:miRNAtotcases
    casefilename = string(miRNAfilenames{i,1});
    casedata = readtable(casefilename);
    casemildata=casedata(:,2);
    casemildata=table2array(casemildata);
    miRNAcasedatamatrix(:,i)=casemildata;
end

mRNAfilenames=mRNAinitialtabledata(2:end,8);
mRNAtotcases=size(mRNAfilenames);
mRNAtotcases=mRNAtotcases(1);
mRNAfirstfilename = string(mRNAfilenames{1,1});
mRNAfirstdata = readtable (mRNAfirstfilename);
mRNAID=mRNAfirstdata(:, 1);
totmRNA=size(mRNAID);
totmRNA=totmRNA(1);


mRNAcasedatamatrix=zeros(totmRNA,mRNAtotcases);

listFrancesca=ls('GoodData Francesca TCGA-KIRC Data'); 

for i=1:mRNAtotcases
    casefilenameFrancesca = string(mRNAfilenames{i,1});
    casedataFrancesca = readtable(casefilenameFrancesca);
    casemildataFrancesca=casedataFrancesca(:,2);
    casemildataFrancesca=table2array(casemildataFrancesca);
    mRNAcasedatamatrix(:,i)=casemildataFrancesca;
end

differentiatedmiRNAsbrian = readtable('DifferentiallyExpressedmiRNAsBrian.txt');

totdifferentiatedmiRNA=size(differentiatedmiRNAsbrian);
totdifferentiatedmiRNA=totdifferentiatedmiRNA(1);

miRNAcasedatamatrixdifferentiated=zeros(totdifferentiatedmiRNA,miRNAtotcases); 
miRNA_ID = table2cell(casedata);
for k2=1:totdifferentiatedmiRNA
    miRNAname2=differentiatedmiRNAsbrian(k2,1);
    C2=table2cell(miRNAname2);
    for kk2=1:totmiRNA

        O= strcmpi(C2 ,miRNA_ID(kk2, 1));
        if O==1
            miRNAcasedatamatrixdifferentiated(k2,:)= miRNAcasedatamatrix(kk2,:); 
        end
    end
end
% miRNAinitialtable2= miRNAinitialtable(:, ismember(miRNA_samples_names, shared_sample_names));
% mRNAinitialtabledata2 = mRNAinitialtabledata(:, ismember(miRNA_samples_names, shared_sample_names));
%  
 [~, locs] = ismember(mRNAID2, EnsemblData.GeneStableIDVersion);
 
 matchingmRNA = EnsemblData.GeneName(locs(locs ~= 0));
 mRNAcasedatamatrix2 = mRNAcasedatamatrix(locs ~= 0, :);
%  
%  
%  [~,b]=ismember(mRNA_sample_names, shared_sample_names)
%  [sorted_varb, sorting_indicesb] = sort(b);
% mRNAinitialtabledata2 = mRNAinitialtabledata(sorted_varb, :);
% % % % % 

%%
% %% 3. Iterate over all of the miRNA-mRNA pairs within the EnsemblData table.
% 
number_of_pairs = 111599
% 
% % set up vectors to hold the r and p values for each pair
r_values = zeros(111599, 1);
p_values = zeros(111599, 1);
% miRNAID2 = miRNA_ID(:, 1);
% mRNAData= EnsemblData2(:, 2);
listmRNA = differentiatedmiRNAsbrian.de_mirnas;
% 
 for ii = 1:number_of_pairs
     % get the names of the miRNA and mRNA in the iith pair
     current_miRNA = mRNATargetbymiRNADiff.Var4(ii);
     current_mRNA = mRNATargetbymiRNADiff.GeneratedAt_2021_01_1613_41_18(ii);
     
     % use the names to get the expression values of the miRNA and mRNA in
     % the iith pair
      expression_of_current_miRNA = miRNAcasedatamatrixdifferentiated(ismember(listmRNA, current_miRNA), :); 
      expression_of_current_mRNA = mRNAcasedatamatrix2(ismember(matchingmRNA, current_mRNA), :); 

     
     % since it's possible that the mRNA in a given miRDIP pair
     % might not be included in the mRNA data from TCGA, or even that it may
     % be included twice  
     % the failsafes: 
     
     if min(size(expression_of_current_mRNA, 1)) == 0
         r_values(ii) = NaN;
         p_values(ii) = NaN;
     else
         if min(size(expression_of_current_mRNA, 1)) > 1
             expression_of_current_mRNA = mean(expression_of_current_mRNA);
         end   
     
         % make your call to corrcoef, and store the output in r_values and
         % p_values. r and p will be returned as 2x2 matrices, (not
         % scalars).
         [r, p] = corrcoef(expression_of_current_miRNA, expression_of_current_mRNA);

         r_values(ii) = r(1, 2);
         p_values(ii) = p(1, 2);
    
      end
 end
 
 mRNATargetbymiRNADiff_r_values = addvars(mRNATargetbymiRNADiff,r_values, 'After', 'Var8');
 mRNATargetbymiRNADiff_rp_values= addvars(mRNATargetbymiRNADiff_r_values, p_values, 'After', 'r_values');
 mRNATargetbymiRNADiff_rp_values2=removevars(mRNATargetbymiRNADiff_rp_values, {'Var3'});
mRNATargetbymiRNADiffNew = rmmissing(mRNATargetbymiRNADiff_rp_values2);
% % removes rows of the table with NaN values in r and p

% % add r_values and p_values as columns to the mRNATargetbymiRNADiff table

% %% 4. BH Correction
pvalues=mRNATargetbymiRNADiffNew.p_values;
fdr = mafdr(pvalues);
fdr2=array2table(fdr);
mRNATargetbymiRNADiffNew2 = [mRNATargetbymiRNADiffNew fdr2];

% mRNATargetbymiRNADiffNew3 = table2array(mRNATargetbymiRNADiffNew2)
mRNATargetbymiRNADiffNew3 = mRNATargetbymiRNADiffNew2(mRNATargetbymiRNADiffNew2.fdr < 0.05, :);

%%
mRNATargetbymiRNADiffNewFinal = mRNATargetbymiRNADiffNew3(mRNATargetbymiRNADiffNew3.r_values < 0, :);
% % perform your Benjamini-Hochberg correction, and get rid of all of the
% % rows of the table that don't pass. 
mRNATargetbymiRNADiffNewFinal.Properties.VariableNames = {'GeneSymbol' 'CaseID' 'miRNA' 'mirDIPscore' 'unknown1' 'unknown2' 'unknown3' 'r_values' ' p_values' 'false_discovery_rate'}
mRNATargetbymiRNADiffNewFinal2= removevars(mRNATargetbymiRNADiffNewFinal,{'CaseID'});
mRNATargetbymiRNADiffNewFinal3= removevars(mRNATargetbymiRNADiffNewFinal2,{'unknown1'});
mRNATargetbymiRNADiffNewFinal4= removevars(mRNATargetbymiRNADiffNewFinal3,{'unknown2'});
mRNATargetbymiRNADiffNewFinal5= removevars(mRNATargetbymiRNADiffNewFinal4,{'unknown3'});

writetable(mRNATargetbymiRNADiffNewFinal5);


