clear
options = optimset('linprog');
options.Display='off';

files_irr = dir('..\Results_irreversibility_considered\*any_kinetic.mat');
files_obj = dir('..\Results_objective\*any_kinetic.mat');

% if e_j=[0 0 1 0 0] being 1 at position j and zero otherwise equals Y'x
% has a solution, j is a strictly stoichiometric BC

% if BC under optimal biomass but not under irreversibility considered
% scenario, the BC is type-II stoichiometric 

% since we have no blocked reactions we cannot check for type-I
% stoichiometric BCs 

% thus all other BCs are stoichiometric
class = cell(length(files_irr),1);

for i=1:length(files_irr)
    disp(i)
    clearvars -except files_irr files_obj i class options
    
    % load data
    R1=load(strcat(files_irr(i).folder,'\',files_irr(i).name));
    R2=load(strcat(files_obj(i).folder,'\',files_obj(i).name));
    
    class{i}=cell(length(R2.B_o{1})+length(R2.B_out{1}),1);
    
    B=[R2.B_o{1} R2.B_out{1}];
    
    for c=1:length(B) % for each balanced complex
        
        % check for strictly stoichiometric BC
        e = zeros(size(R2.MODEL_o{1}.Y,2),1);
        e(B(c)) = 1;
        
        Sol=linprog(zeros(size(R2.MODEL_o{1}.Y',2),1),[],[],R2.MODEL_o{1}.Y',e,[],[],[],options);
        if ~isempty(Sol)
            class{i}(c) = {'stricktly stoichiometric'};
        
        % check for type-II nonstoichiometric BC
        elseif isempty(intersect([R1.B_r{1} R1.B_out{1}],B(c)))
            class{i}(c) = {'type-II nonstoichiometric'};
            
        else
            class{i}(c) = {'stoichiometric'};
        end          
    end     
end


for i=1:length(class)
    stricktly_stoich(i,1) = sum(strcmp(class{i},'stricktly stoichiometric'))/length(class{i});
    stoich(i,1) = sum(strcmp(class{i},'stoichiometric'))/length(class{i});
    typeII(i,1) = sum(strcmp(class{i},'type-II nonstoichiometric'))/length(class{i});
    label(i,1)={files_irr(i).name(1:end-16)};
end

label=strrep(label,'_',' ')

label={'\itA. niger iMA871';
    '\itA. thaliana AraCore';
    '\itC. reinhardtii iCre1355';
    '\itE. coli iJO1366';
    '\itM. acetivorans iMB745';
    '\itM. barkeri iAF692';
    '\itM. musculus';
    '\itM. tuberculosis iNJ661m';
    '\itN. pharaonis';
    '\itP. putida iJN746';
    '\itT. maritima iLJ478';
    '\itS. cerevisiae Yeast8'};
    
bar([stricktly_stoich stoich typeII]*100,'stacked')
set(gca,'XTick',1:12,'XTickLabel',label,'XTickLabelRotation',45)
ylabel('Distribution of balanced complex fractorizations (%)')
legend('strictly stoichiometric','stoichiometric','type-II nonstoichiometric',...
    'Location','NorthEastOutside')
legend boxoff