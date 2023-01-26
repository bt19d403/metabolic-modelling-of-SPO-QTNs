clear
initCobraToolbox(false)
changeCobraSolver('gurobi','all')
model = readCbModel('yeastGEM.mat')


model = changeRxnBounds(model,'r_1106',-1000,'l'); % acetate[e] to acetate [c]
model = changeRxnBounds(model,'r_1106',0,'u'); % acetate[e] to acetate [c]
model = changeRxnBounds(model,'r_1634',-1000,'l'); % spitting out/in acetate[e]
model = changeRxnBounds(model,'r_1634',0,'u'); % spitting out/in acetate[e]

model = changeRxnBounds(model,'r_1992',-1000,'l'); % oxygen uptake
model = changeRxnBounds(model,'r_1979',-1000,'l'); % oxygen uptake

model = changeRxnBounds(model,'r_1714',0,'l'); % glucose uptake
model = changeRxnBounds(model,'r_1115',0,'l'); % ammonia


% % for biomass growth
% model.c(:)=0;
% model.c(2906) = 1;
% f = optimizeCbModel(model);
% bio_flux = f.f;
% % for biomass maintanence
% model.c(:)=0;
% model.c(3414) = 1;
% f = optimizeCbModel(model);
% bio_main_flux = f.f;
% % for glut_exc
% model.c(:)=0;
% model.c(2737) = 1;
% f = optimizeCbModel(model);
% glu_exc_flux = f.f;
% % for nuc
% model.c(:)=0;
% model.c(361) = 1;
% f = optimizeCbModel(model);
% nuc_flux = f.f;
% 

%% get IMAT scores using localginni
geneExpData = readtable('tpm_counts_Average.csv');
expData.value = geneExpData(:,2:17);
expData.value = table2array(geneExpData(:,2:17));
expData.genes = table2cell(geneExpData(:,1));
geneExpData.Properties.VariableNames
context = geneExpData.Properties.VariableNames;
context(:,1) = [];
expData.context = context
%model extraction method
MeM = 'iMAT'
contexts = {'C1','C2','C3','C4','C5','C6','C7','C8','C9','C10','C11','C12','C13','C14','C15','C16'};
ut = 90;
lt = 10;
ThS = 1; % impliying at gene level
%Tolerance level above which reactions are considered as expressed
tol = 1e-4;
filename = './';
biomass_id = find(strcmp(model.rxns,'r_2111'));
nucleotide_synthesis_id = find(strcmp(model.rxns,'r_0466'));
NGAM_id = find(strcmp(model.rxns,'r_4046'));
% model.c(:)=0;
% model.c(biomass_id)=1;
% model.ub(find(model.c))=0.02*3.9943; % for 2 percent of max growth
% result = optimizeCbModel(model,'max','one');
% rxn = find(result.x ~= 0);
% coreRxn=[biomass_id , nucleotide_synthesis_id, NGAM_id];
coreRxn=[biomass_id , nucleotide_synthesis_id, NGAM_id]
[Models,RxnImp] = buildContextmodels(expData,model,MeM,contexts,ut,lt,ThS,coreRxn,filename,tol);
% 
% 
%

 for i=1:16
    temp_var = removeUnusedGenes(Models{i});
    Models{i} = temp_var;
end

C1_model = Models{1}
C2_model = Models{2}
C3_model = Models{3}
C4_model = Models{4}
C5_model = Models{5}
C6_model = Models{6}
C7_model = Models{7}
C8_model = Models{8}
C9_model = Models{9}
C10_model = Models{10}
C11_model = Models{11}
C12_model = Models{12}
C13_model = Models{13}
C14_model = Models{14}
C15_model = Models{15}
C16_model = Models{16}

%% 
for i = 1:16
[Nucleotide_flux(i) ,num_rxns(i), num_genes(i), num_metabolites(i)] = find_nucleotide_flux(Models{i})

end

%% 
for i = 1:16
[Ngam(i)] = find_ngam_flux(Models{i});

end

%% 
rxns_set1 = C1_model.rxns;
rxns_set2 = C2_model.rxns;
rxns_set3 = C3_model.rxns;
rxns_set4 = C4_model.rxns;
rxns_set5 = C5_model.rxns;
rxns_set6= C6_model.rxns;
rxns_set7 = C7_model.rxns;
rxns_set8= C8_model.rxns;
rxns_set9= C9_model.rxns;
rxns_set10 = C10_model.rxns;
rxns_set11 = C11_model.rxns;
rxns_set12 = C12_model.rxns;
rxns_set13 = C13_model.rxns;
rxns_set14 = C14_model.rxns;
rxns_set15 = C15_model.rxns;
rxns_set16 = C16_model.rxns;


% unique reactions including wildtype
unique_rxns_1 = find(~ismember(rxns_set1,rxns_set16));
unique_rxns_12 =find(~ismember(rxns_set12,rxns_set16));
unique_rxns_9 =find(~ismember(rxns_set9,rxns_set16));
unique_rxns_15 =find(~ismember(rxns_set15,rxns_set16));
unique_rxns_2 =find(~ismember(rxns_set2,rxns_set16));
unique_rxns_3 = find(~ismember(rxns_set3,rxns_set16)); 
unique_rxns_4 = find(~ismember(rxns_set4,rxns_set16)); 
unique_rxns_5 = find(~ismember(rxns_set5,rxns_set16)); 
unique_rxns_6 = find(~ismember(rxns_set6,rxns_set16)); 
unique_rxns_7 = find(~ismember(rxns_set7,rxns_set16));
unique_rxns_8 = find(~ismember(rxns_set8,rxns_set16)); 
unique_rxns_10= find(~ismember(rxns_set10,rxns_set16)); 
unique_rxns_11 = find(~ismember(rxns_set11,rxns_set16)); 
unique_rxns_13 = find(~ismember(rxns_set13,rxns_set16)); 
unique_rxns_14 = find(~ismember(rxns_set14,rxns_set16)); 


%% FLux enrichment for unique reactions
resultCell_C1_model_unique = FEA(C1_model, unique_rxns_1, 'subSystems');
resultCell_C2_model_unique = FEA(C2_model, unique_rxns_2, 'subSystems');
resultCell_C3_model_unique = FEA(C3_model, unique_rxns_3, 'subSystems');
resultCell_C4_model_unique = FEA(C4_model, unique_rxns_4, 'subSystems');
resultCell_C5_model_unique = FEA(C5_model, unique_rxns_5, 'subSystems');
resultCell_C6_model_unique = FEA(C6_model, unique_rxns_6, 'subSystems');
resultCell_C7_model_unique = FEA(C7_model, unique_rxns_7, 'subSystems');
resultCell_C8_model_unique = FEA(C8_model, unique_rxns_8, 'subSystems');
resultCell_C9_model_unique = FEA(C9_model, unique_rxns_9, 'subSystems');
resultCell_C10_model_unique = FEA(C10_model, unique_rxns_10, 'subSystems');
resultCell_C11_model_unique = FEA(C11_model, unique_rxns_11, 'subSystems');
resultCell_C12_model_unique = FEA(C12_model, unique_rxns_12, 'subSystems');
resultCell_C13_model_unique = FEA(C13_model, unique_rxns_13, 'subSystems');
resultCell_C14_model_unique = FEA(C14_model, unique_rxns_14, 'subSystems');
resultCell_C15_model_unique = FEA(C15_model, unique_rxns_15, 'subSystems');


%% 

C1_model = readCbModel('C1_model_new_lg_imat.mat');
C2_model = readCbModel('C2_model_new_lg_imat.mat')
C3_model = readCbModel('C3_model_new_lg_imat.mat')
C4_model = readCbModel('C4_model_new_lg_imat.mat')
C5_model = readCbModel('C5_model_new_lg_imat.mat')
C6_model = readCbModel('C6_model_new_lg_imat.mat')
C7_model = readCbModel('C7_model_new_lg_imat.mat')
C8_model = readCbModel('C8_model_new_lg_imat.mat')
C9_model = readCbModel('C9_model_new_lg_imat.mat')
C10_model = readCbModel('C10_model_new_lg_imat.mat')
C11_model = readCbModel('C11_model_new_lg_imat.mat')
C12_model = readCbModel('C12_model_new_lg_imat.mat')
C13_model = readCbModel('C13_model_new_lg_imat.mat')
C14_model = readCbModel('C14_model_new_lg_imat.mat')
C16_model = readCbModel('C16_model_new_lg_imat.mat')
C15_model = readCbModel('C15_model_new_lg_imat.mat')

%%% chnage objective to ngam 



% ngam objective

C1_model_obj = changeObjective(C1_model,'r_4046');
C2_model_obj = changeObjective(C2_model,'r_4046');
C3_model_obj = changeObjective(C3_model,'r_4046');
C4_model_obj = changeObjective(C4_model,'r_4046');
C5_model_obj = changeObjective(C5_model,'r_4046');
C6_model_obj = changeObjective(C6_model,'r_4046');
C7_model_obj = changeObjective(C7_model,'r_4046');
C8_model_obj = changeObjective(C8_model,'r_4046');
C9_model_obj = changeObjective(C9_model,'r_4046');
C10_model_obj = changeObjective(C10_model,'r_4046');
C11_model_obj = changeObjective(C11_model,'r_4046');
C12_model_obj = changeObjective(C12_model,'r_4046');
C13_model_obj = changeObjective(C13_model,'r_4046');
C14_model_obj = changeObjective(C14_model,'r_4046');
C15_model_obj = changeObjective(C15_model,'r_4046');
C16_model_obj = changeObjective(C16_model,'r_4046');

C1_model_obj.ub(find(C1_model_obj.c)) = 1000;
C2_model_obj.ub(find(C2_model_obj.c)) = 1000;
C3_model_obj.ub(find(C3_model_obj.c)) = 1000;
C4_model_obj.ub(find(C4_model_obj.c)) = 1000;
C5_model_obj.ub(find(C5_model_obj.c)) = 1000;
C6_model_obj.ub(find(C6_model_obj.c)) = 1000;
C7_model_obj.ub(find(C7_model_obj.c)) = 1000;
C8_model_obj.ub(find(C8_model_obj.c)) = 1000;
C9_model_obj.ub(find(C9_model_obj.c)) = 1000;
C10_model_obj.ub(find(C10_model_obj.c)) = 1000;
C11_model_obj.ub(find(C1_model_obj.c)) = 1000;
C1_model_obj.ub(find(C11_model_obj.c)) = 1000;
C12_model_obj.ub(find(C12_model_obj.c)) = 1000;
C13_model_obj.ub(find(C13_model_obj.c)) = 1000;
C14_model_obj.ub(find(C14_model_obj.c)) = 1000;
C15_model_obj.ub(find(C15_model_obj.c)) = 1000;
C16_model_obj.ub(find(C16_model_obj.c)) = 1000;
%flux variability
%% 

[minFlux1, maxFlux1] = fluxVariability(C1_model_obj,2);
[minFlux16, maxFlux16] = fluxVariability(C16_model_obj,2);
[minFlux2, maxFlux2] = fluxVariability(C2_model_obj,2);
[minFlux15, maxFlux15] = fluxVariability(C15_model_obj,2);
[minFlux12, maxFlux12] = fluxVariability(C12_model_obj,2);
[minFlux9, maxFlux9] = fluxVariability(C9_model_obj,2);
[minFlux3, maxFlux3] = fluxVariability(C3_model_obj,2);
[minFlux4, maxFlux4] = fluxVariability(C4_model_obj,2);
[minFlux5, maxFlux5] = fluxVariability(C5_model_obj,2);
[minFlux6, maxFlux6] = fluxVariability(C6_model_obj,2);
[minFlux7, maxFlux7] = fluxVariability(C7_model_obj,2);
[minFlux8, maxFlux8] = fluxVariability(C8_model_obj,2);
[minFlux9, maxFlux9] = fluxVariability(C9_model_obj,2);
[minFlux10, maxFlux10] = fluxVariability(C10_model_obj,2);
[minFlux11, maxFlux11] = fluxVariability(C11_model_obj,2);
[minFlux13, maxFlux13] = fluxVariability(C13_model_obj,2);
[minFlux14, maxFlux14] = fluxVariability(C14_model_obj,2);

%% identify common reactions between all models

Common_all_models = readtable('common_rxns_c16_new.csv');
resultCell_common = FEA(C1_model, Common_all_models.('Reaction_num1'), 'subSystems');

%%  identify common reactions between all SNP  models

Common_all_SNPmodels = readtable('common_rxns_SNP_models.csv');
resultCell_common_SNP = FEA(C1_model, Common_all_SNPmodels.('Reaction_num1'), 'subSystems');
%% 
initCobraToolbox(false);
changeCobraSolver('gurobi','all')
C1_model = readCbModel('C1_model_new_lg_imat.mat');
C2_model = readCbModel('C2_model_new_lg_imat.mat')
C3_model = readCbModel('C3_model_new_lg_imat.mat')
C4_model = readCbModel('C4_model_new_lg_imat.mat')
C5_model = readCbModel('C5_model_new_lg_imat.mat')
C6_model = readCbModel('C6_model_new_lg_imat.mat')
C7_model = readCbModel('C7_model_new_lg_imat.mat')
C8_model = readCbModel('C8_model_new_lg_imat.mat')
C9_model = readCbModel('C9_model_new_lg_imat.mat')
C10_model = readCbModel('C10_model_new_lg_imat.mat')
C11_model = readCbModel('C11_model_new_lg_imat.mat')
C12_model = readCbModel('C12_model_new_lg_imat.mat')
C13_model = readCbModel('C13_model_new_lg_imat.mat')
C14_model = readCbModel('C14_model_new_lg_imat.mat')
C16_model = readCbModel('C16_model_new_lg_imat.mat')
C15_model = readCbModel('C15_model_new_lg_imat.mat')
%% 
up_c1 = readtable('up_c1_16_ngam.csv');
up_c2 = readtable('up_c2_16_ngam.csv');
up_c3 = readtable('up_c3_16_ngam.csv');
up_c4 = readtable('up_c4_16_ngam.csv');
up_c5 = readtable('up_c5_16_ngam.csv');
up_c6 = readtable('up_c6_16_ngam.csv');
up_c7 = readtable('up_c7_16_ngam.csv');
up_c8 = readtable('up_c8_16_ngam.csv');
up_c9 = readtable('up_c9_16_ngam.csv');
up_c10 = readtable('up_c10_16_ngam.csv');
up_c11= readtable('up_c11_16_ngam.csv');
up_c12= readtable('up_c12_16_ngam.csv');
up_c13= readtable('up_c13_16_ngam.csv');
up_c14= readtable('up_c14_16_ngam.csv');
up_c15= readtable('up_c15_16_ngam.csv');


%% 
resultCell_C1_model_up = FEA(C16_model, up_c1.('Reaction_num16'), 'subSystems');
resultCell_C2_model_up = FEA(C16_model, up_c2.('Reaction_num16'), 'subSystems');
resultCell_C3_model_up = FEA(C16_model, up_c3.('Reaction_num16'), 'subSystems');
resultCell_C4_model_up = FEA(C16_model, up_c4.('Reaction_num16'), 'subSystems');
resultCell_C5_model_up = FEA(C16_model, up_c5.('Reaction_num16'), 'subSystems');
resultCell_C6_model_up = FEA(C16_model, up_c6.('Reaction_num16'), 'subSystems');
resultCell_C7_model_up = FEA(C16_model, up_c7.('Reaction_num16'), 'subSystems');
resultCell_C8_model_up = FEA(C16_model, up_c8.('Reaction_num16'), 'subSystems');
resultCell_C9_model_up = FEA(C16_model, up_c9.('Reaction_num16'), 'subSystems');
resultCell_C10_model_up = FEA(C16_model, up_c10.('Reaction_num16'), 'subSystems');
resultCell_C11_model_up = FEA(C16_model, up_c11.('Reaction_num16'), 'subSystems');
resultCell_C12_model_up = FEA(C16_model, up_c12.('Reaction_num16'), 'subSystems');
resultCell_C13_model_up = FEA(C16_model, up_c13.('Reaction_num16'), 'subSystems');
resultCell_C14_model_up = FEA(C16_model, up_c14.('Reaction_num16'), 'subSystems');
resultCell_C15_model_up = FEA(C16_model, up_c15.('Reaction_num16'), 'subSystems');

%% 
down_c1 = readtable('down_c1_16_ngam.csv');
down_c2 = readtable('down_c2_16_ngam.csv');
down_c3 = readtable('down_c3_16_ngam.csv');
down_c4 = readtable('down_c4_16_ngam.csv');
down_c5 = readtable('down_c5_16_ngam.csv');
down_c6 = readtable('down_c6_16_ngam.csv');
down_c7 = readtable('down_c7_16_ngam.csv');
down_c8 = readtable('down_c8_16_ngam.csv');
down_c9 = readtable('down_c9_16_ngam.csv');
down_c10 = readtable('down_c10_16_ngam.csv');
down_c11= readtable('down_c11_16_ngam.csv');
down_c12= readtable('down_c12_16_ngam.csv');
down_c13= readtable('down_c13_16_ngam.csv');
down_c14= readtable('down_c14_16_ngam.csv');
down_c15= readtable('down_c15_16_ngam.csv');

%% 
resultCell_C1_model_down = FEA(C16_model, down_c1.('Reaction_num16'), 'subSystems');
resultCell_C2_model_down = FEA(C16_model, down_c2.('Reaction_num16'), 'subSystems');
resultCell_C3_model_down = FEA(C16_model, down_c3.('Reaction_num16'), 'subSystems');
resultCell_C4_model_down = FEA(C16_model, down_c4.('Reaction_num16'), 'subSystems');
resultCell_C5_model_down = FEA(C16_model, down_c5.('Reaction_num16'), 'subSystems');
resultCell_C6_model_down = FEA(C16_model, down_c6.('Reaction_num16'), 'subSystems');
resultCell_C7_model_down = FEA(C16_model, down_c7.('Reaction_num16'), 'subSystems');
resultCell_C8_model_down = FEA(C16_model, down_c8.('Reaction_num16'), 'subSystems');
resultCell_C9_model_down = FEA(C16_model, down_c9.('Reaction_num16'), 'subSystems');
resultCell_C10_model_down = FEA(C16_model, down_c10.('Reaction_num16'), 'subSystems');
resultCell_C11_model_down = FEA(C16_model, down_c11.('Reaction_num16'), 'subSystems');
resultCell_C12_model_down = FEA(C16_model, down_c12.('Reaction_num16'), 'subSystems');
resultCell_C13_model_down = FEA(C16_model, down_c13.('Reaction_num16'), 'subSystems');
resultCell_C14_model_down = FEA(C16_model, down_c14.('Reaction_num16'), 'subSystems');
resultCell_C15_model_down = FEA(C16_model, down_c15.('Reaction_num16'), 'subSystems');
%% 

reaction_name_C1_model_up =  C16_model.rxnNames(up_c1.('Reaction_num16'))
reaction_name_C2_model_up = C16_model.rxnNames(up_c2.('Reaction_num16'))
reaction_name_C3_model_up = C16_model.rxnNames(up_c3.('Reaction_num16'))
reaction_name_C4_model_up = C16_model.rxnNames(up_c4.('Reaction_num16'))
reaction_name_C5_model_up = C16_model.rxnNames(up_c5.('Reaction_num16'))
reaction_name_C6_model_up = C16_model.rxnNames(up_c6.('Reaction_num16'))
reaction_name_C7_model_up = C16_model.rxnNames(up_c7.('Reaction_num16'))
reaction_name_C8_model_up = C16_model.rxnNames(up_c8.('Reaction_num16'))
reaction_name_C9_model_up = C16_model.rxnNames(up_c9.('Reaction_num16'))
reaction_name_C10_model_up = C16_model.rxnNames(up_c10.('Reaction_num16'))
reaction_name_C11_model_up = C16_model.rxnNames(up_c11.('Reaction_num16'))
reaction_name_C12_model_up = C16_model.rxnNames(up_c12.('Reaction_num16'))
reaction_name_C13_model_up =C16_model.rxnNames(up_c13.('Reaction_num16'))
reaction_name_C14_model_up = C16_model.rxnNames(up_c14.('Reaction_num16'))
reaction_name_C15_model_up = C16_model.rxnNames(up_c15.('Reaction_num16'))

%% 
reaction_name_C1_model_down =  C16_model.rxnNames(down_c1.('Reaction_num16'))
reaction_name_C2_model_down = C16_model.rxnNames(down_c2.('Reaction_num16'))
reaction_name_C3_model_down = C16_model.rxnNames(down_c3.('Reaction_num16'))
reaction_name_C4_model_down = C16_model.rxnNames(down_c4.('Reaction_num16'))
reaction_name_C5_model_down = C16_model.rxnNames(down_c5.('Reaction_num16'))
reaction_name_C6_model_down = C16_model.rxnNames(down_c6.('Reaction_num16'))
reaction_name_C7_model_down = C16_model.rxnNames(down_c7.('Reaction_num16'))
reaction_name_C8_model_down = C16_model.rxnNames(down_c8.('Reaction_num16'))
reaction_name_C9_model_down = C16_model.rxnNames(down_c9.('Reaction_num16'))
reaction_name_C10_model_down = C16_model.rxnNames(down_c10.('Reaction_num16'))
reaction_name_C11_model_down = C16_model.rxnNames(down_c11.('Reaction_num16'))
reaction_name_C12_model_down = C16_model.rxnNames(down_c12.('Reaction_num16'))
reaction_name_C13_model_down =C16_model.rxnNames(down_c13.('Reaction_num16'))
reaction_name_C14_model_down = C16_model.rxnNames(down_c14.('Reaction_num16'))
reaction_name_C15_model_down = C16_model.rxnNames(down_c15.('Reaction_num16'))

