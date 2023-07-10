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




