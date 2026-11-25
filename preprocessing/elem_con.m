function ub = elem_con(model,metCarb,metNitro,bnds,varID)

%% extract lower and upper bnds from bnds variable
if iscell(bnds)
    k = 0;
    [num_vars,~] = size(bnds);
    for i=1:num_vars
        k = k+1;
        varMin(k,1) = bnds{i,2};
    end
end

%% find corresponding metIDs from varIDs
% first remove all bnd rxn that are not true drains (e.g. biomass, IgG formation)
isDrain = find(checkDrainRxns(model));
noExRxns = find(~ismember(varID,isDrain));
% adjust the previously retrieved variables
varID(noExRxns) = [];
varMin(noExRxns) = [];
% get all drains without the set bnds
isDrainNew = isDrain(~ismember(isDrain,varID));
% retrieve the actual metabolite IDs for carbon calculations
metID_tmp = find(model.S(:,varID));
[metID,~] = ind2sub(size(model.S),metID_tmp); % index to coordinates

%% multiply the lower and upper bnds witht their respective number of carbon atoms
varMinCarb = varMin.*metCarb(metID);
varMinNitro = varMin.*metNitro(metID);
total_carbon_cons = abs(sum(varMinCarb(varMinCarb<0)));
total_nitrogen_cons = abs(sum(varMinNitro(varMinNitro<0)));

    
%% Constrain biomass production
% c = 0;
% for i = 1:num_rxns    
%     if ismember(i,isDrain)        
%         % Carbon constrain
%         subIDs = find(model.S(:,i) < 0);
%         prodIDs = find(model.S(:,i) > 0);
%         subCarbon = sum(abs(full(model.S(subIDs,i))).*metCarb(subIDs,1));
%         prodCarbon = sum(abs(full(model.S(prodIDs,i))).*metCarb(prodIDs,1));
%         if subCarbon == 0
%             carbon_constraint = abs(total_carbon_cons/prodCarbon);
%         else
%             carbon_constraint = abs(total_carbon_cons/subCarbon);
%         end
%         
%         % Nitrogen Constrain
%         subNitrogen = sum(abs(full(model.S(subIDs,i))).*metNitro(subIDs,1));
%         prodNitrogen = sum(abs(full(model.S(prodIDs,i))).*metNitro(prodIDs,1));
%         if subNitrogen == 0
%             nitrogen_constraint = abs(total_nitrogen_cons/prodNitrogen);
%         else
%             nitrogen_constraint = abs(total_nitrogen_cons/subNitrogen);
%         end
%         
%         if model.rev(i,1) == 1
%             constraint = min(carbon_constraint,nitrogen_constraint);
%             c = c + 1;
%             ub(c) = constraint;
%         end
%     end
% end

% Biomass
subIDs = find(model.S(:,5190) < 0);
subCarbon = sum(abs(full(model.S(subIDs,5190))).*metCarb(subIDs,1));
subNitrogen = sum(abs(full(model.S(subIDs,5190))).*metNitro(subIDs,1));

carbon_constraint = abs(total_carbon_cons/subCarbon);
nitrogen_constraint = abs(total_nitrogen_cons/subNitrogen);
constraint = min(carbon_constraint,nitrogen_constraint);
ub(1) = constraint;

% IgG
subIDs = find(model.S(:,5195) < 0);
subCarbon = sum(abs(full(model.S(subIDs,5195))).*metCarb(subIDs,1));
subNitrogen = sum(abs(full(model.S(subIDs,5195))).*metNitro(subIDs,1));
constraint = min(carbon_constraint,nitrogen_constraint);

carbon_constraint = abs(total_carbon_cons/subCarbon);
nitrogen_constraint = abs(total_nitrogen_cons/subNitrogen);
constraint = min(carbon_constraint,nitrogen_constraint);

ub(2) = constraint;

subIDs = find(model.S(:,5196) < 0);
subCarbon = sum(abs(full(model.S(subIDs,5196))).*metCarb(subIDs,1));
subNitrogen = sum(abs(full(model.S(subIDs,5196))).*metNitro(subIDs,1));
constraint = min(carbon_constraint,nitrogen_constraint);

carbon_constraint = abs(total_carbon_cons/subCarbon);
nitrogen_constraint = abs(total_nitrogen_cons/subNitrogen);
constraint = min(carbon_constraint,nitrogen_constraint);

ub(3) = constraint;

subIDs = find(model.S(:,5197) < 0);
subCarbon = sum(abs(full(model.S(subIDs,5197))).*metCarb(subIDs,1));
subNitrogen = sum(abs(full(model.S(subIDs,5197))).*metNitro(subIDs,1));
constraint = min(carbon_constraint,nitrogen_constraint);

carbon_constraint = abs(total_carbon_cons/subCarbon);
nitrogen_constraint = abs(total_nitrogen_cons/subNitrogen);
constraint = min(carbon_constraint,nitrogen_constraint);

ub(4) = constraint;

end