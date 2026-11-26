function [carbonConst] = carbonConstraints(model,bnds,condition,carbonCount,cofPairs,minmax,runFinalMinMax,limEX,totalCarbon,skipRxns)
[~,num_rxns] = size(model.S);

%% This script is designed to constrain the solution space based on the 
%   number of carbon atoms made available to the model by the exchange
%   reactions. This is done in order to avoid physiologically meaningless
%   futile cycles.
%
%   CREATED 2017-09-30 by Maximilian Lularevic
%   Edited by Athanasios Antonakoudis 10/12/2018
%
% INPUT
%   model       =   cobra model (must contain model.metFormulas)
%
%   bnds        =   bounds or constraints as cell or double. for double 
%                   first row corresponds to rxnIDs, second row to lower 
%                   and the third row to upper bnds
%
% OPTIONAL
%
%   condition   =   string containing information as to which condition is 
%                   run (e.g. lactate producing conditon)
%
%   carbonCount =   carbonCount Structure as received form function 
%                   findMetCarbon.m (needs to contain carbonCount.carbonAtoms field)
%
%   cofParis    =   cofParis structure as received from function findCofactorPairs.m
%
%   minmax      =   FVA for model under experimental constraints. (needs 
%                   to be processed with fixMinMax.m)
%
%   runFinalMM  =   logical true or false
%
%   limEX       =   value to be used to constrain drain rxns. This is
%                   supposed to be a single carbon constraint (e.g. left
%                   over carbon can be set to CO2 rate which is a single
%                   carbon. If flux for CO2 (closing the carbon in excel)
%                   is 0.641 mmol/gDCW*h this would be the rate to fed into
%                   the limEX variable).
%
% OUTPUT
%   carbonConst =   structure containing the following fields
%                       - Title
%                       - Date and Time of job
%                       - constraintsName (given: e.g. 'LacProd')
%                       - model.description
%                       - total carbon consumed
%               *** optional ***
%                       - actual carbon used for constraining (only when
%                         passed into function)
%               *** optional ***
%                       - report
%                           - cell with the lenght of model.rxns which 
%                             states whther reaction was balanced or unbalanced  
%                           - number of balanced rxns
%                           - number of unbalanced rxns
%                           - number of set rxns (aka constraints)
%                           - number of rxn with no flux
%                       - bndAdjusted
%                           - rawData (logical 2 colum matrix)
%                           - number of lower bounds adjusted by cc
%                             algorithm
%                           - number of upper bnds adjusted by cc algorithm
%                       - carbonCount
%                           - dateAndTime
%                           - modelDescription
%                           - metName
%                           - metFormula
%                           - carbon atoms for each metabolite
%                           - table with all the above information compiled
%                       - cofactorPairs (cell with mets and met IDs and
%                         occurance of these mets with corresponding
%                         rxnIDs)
%                       - original minmax (either the one passed on to 
%                         function or herein calculated one)
%                       - newMinMax pased on carbon balance
%                       - solSpaceReduction size of new ssolution space
%                         compared to old one (calculatons based on minmax
%                         range comparison)
%                       - solutionNew (optimizeCbModel output with new bnds)
%                       - EXrxnsAdjusted -> title
%                       - limEX (single carbon limit passed into fucntion
%                         to constrain the Exchange reactions output)
%                       - number of EX bnds adjusted
%                       - solSpaceReductionEX in % compared to original
%                       - EXspaceReduction (how much has the space shrunk
%                         just comparing old EX ranges with new ones
%                       - post New MinMax -> title
%                       - postNewMinMax (new minmax run with runMinMax_GF
%                         and cc bnds)
%                       - solSpacereductionPostCC
%                       - solutionPostNew


%% check if minmax was passed into function and, run minmax if it hasn't
if exist('minmax','var')
    if isempty(minmax)
        minmax = runMinMax_GF(model);
        minmax = fixMinMax(minmax);
    end
else
    max = 1000;
    model.rev=checkReversability(model);
    model = resetModel(model,-max,max,model.rev);
    [lb,ub,model] = Medium_Bounds(model,max); 
    model = set_ub(model,ub);
    model = set_lb(model,lb);
    [lb,ub] = Susie_Bounds(1);
    model = set_ub(model,ub);
    model = set_lb(model,lb);
    minmax = runMinMax_GF(model); % runs the simple MinMax algorithm for the model;
    minmax = fixMinMax(minmax);   % fixes the MinMax, if any lb>ub; usually there aren't any
 end

%% check if condition was passed into function
if exist('condition','var')
    if isempty(condition)
        condition = 'n/a';
    end
else
    condition = 'n/a';
end

%% check if there are any rxns to be skipped
if exist('skipRxns','var')
    if isempty(skipRxns)
        skipRxns = [];
    else
        [s,z] = size(skipRxns);
        if z > s
            skipRxns = skipRxns';
        end
    end
else
    skipRxns = [];
end

%% extract lower and upper bnds from bnds variable
if iscell(bnds)
    k = 0;
    [num_vars,~] = size(bnds);
    for i=1:num_vars
        varName = bnds{i,1};
        var_index = find(ismember(model.rxns,varName));
        if isempty(var_index)
            fprintf('%s not found\n',bnds{i,1});
        else
            k = k+1;
            varID(k,1) = var_index;
            varMin(k,1) = bnds{i,2};
            varMax(k,1) = bnds{i,3};
        end
    end
elseif isnumeric(bnds)
    varID = bnds(:,1);
    varMin = bnds(:,2);
    varMax = bnds(:,3);
end
varID_tmp = varID;

%% find corresponding metIDs from varIDs
% first remove all bnd rxn that are not true drains (e.g. biomass, IgG formation)
isDrain = find(checkDrainRxns(model));
noExRxns = find(~ismember(varID,isDrain));
% adjust the previously retrieved variables
varID(noExRxns) = [];
varMin(noExRxns) = [];
varMax(noExRxns) = [];
% get all drains without the set bnds
isDrainNew = isDrain(~ismember(isDrain,varID));
% retrieve the actual metabolite IDs for carbon calculations
metID_tmp = find(model.S(:,varID));
[metID,~] = ind2sub(size(model.S),metID_tmp); % index to coordinates

%% check if number of carbon atoms per metabolite was passed into function
if exist('carbonCount','var')
    if isempty(carbonCount)
        carbonCount = findMetCarbon(model);
    end
else
    carbonCount = findMetCarbon(model);
end
% extract the carbon vector (length of model.mets)
metCarb = carbonCount.carbonAtoms;

%% multiply the lower and upper bnds witht their respective number of carbon atoms
varMinCarb = varMin.*metCarb(metID);
if exist('totalCarbon','var')
    if isempty(totalCarbon)
        total_carbon_cons = abs(sum(varMinCarb(varMinCarb<0)));
    elseif ~isnumeric(totalCarbon)
        h = warndlg('Warning the totalCarbon input is not numeric. The totalCarbon will be calculated based on the lower bnds of the input variables. If you wish to continue press OK or else kill the code and input a numeric value for totalCarbon.','Problem with input');
        total_carbon_cons = abs(sum(varMinCarb(varMinCarb<0)));
    else
        total_carbon_cons = totalCarbon;
        actual_total_carbon_cons = abs(sum(varMinCarb(varMinCarb<0)));
    end
else
    total_carbon_cons = abs(sum(varMinCarb(varMinCarb<0)));

end

%% check if cofactorPairs have been passed on and set cofactors to 0 in model.S

% check for cofPari variable and/or calculate new
if exist('cofPairs','var')
    if isempty(cofPairs)
        smallMets = findSmallMets(model);
        metPairs = countMetPairs_max(model,smallMets.smallMetIDs,true,'[');
        cofPairs = findCofactorPairs_max(metPairs);
    end
else
    smallMets = findSmallMets(model);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%% THIS NEEDS CHANGE %%%%%%%%%%%%%%%%%%%%%%%%%%
    metPairs = countMetPairs_max(model,smallMets.smallMetIDs,true,'_');
    %%%%%%%%%%%%%%%%%%%%%%%%% THIS NEEDS CHANGE %%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cofPairs = findCofactorPairs_max(metPairs);
end

% initialize variable
cofMetIDs = [];
% retrieve all the metabolite IDs for cofactors
for i = 1:size(cofPairs,1)
    cof_tmp = cell2mat(cofPairs{i,6});
    cofMetIDs = [cofMetIDs;cof_tmp];
end

% set cofMets to 0 in model.S
model_tmp = model;
model_tmp.S(cofMetIDs,:) = 0;

% find all rxns containing coa mets on both sides of the rxn
coa = findCoAs(model_tmp,' ','C');

% find all rxns containing cytochromes, quinones, and hemes
lrgMlc = findLrgMlcs(model,metCarb);

%% calculating new bnds algorithm !!!!!
minmax_new = minmax;
mmAdj = zeros(num_rxns,2);
for i = 1:num_rxns


    % First, check if the rxn is an important one
    if (i == find(ismember(model.rxns,'biomass_cho_producing'))) % Biomass non producing
        % Get the coefficients for the biomass, not using the model_temp
        % because we want the coefpairs
        subIDs = find(model.S(:,i) < 0);
        % Calculating the carbon same as b4
        subCarbon = sum(abs(full(model.S(subIDs,i))).*metCarb(subIDs,1));
        bericht{i,1} ='Biomass nonproducing';
        constraint = abs(total_carbon_cons/subCarbon);
        if abs(minmax(i,2)) > constraint
            minmax_new(i,2) = constraint;
            mmAdj(i,2) = 1;
        end
        continue
    elseif (i == find(ismember(model.rxns,'biomass_cho'))) %Biomass producing
        % Get the coefficients for the biomass, not using the model_temp
        % because we want the coefpairs
        subIDs = find(model.S(:,i) < 0);
        % Calculating the carbon same as b4
        subCarbon = sum(abs(full(model.S(subIDs,i))).*metCarb(subIDs,1));
        bericht{i,1} ='Biomass producing';
            constraint = abs(total_carbon_cons/subCarbon);
            if abs(minmax(i,2)) > constraint
                minmax_new(i,2) = constraint;
                mmAdj(i,2) = 1;
            end
            continue
    elseif (i == find(ismember(model.rxns,'igg_hc'))) %HC Synthesis
        % Get the coefficients for the biomass, not using the model_temp
        % because we want the coefpairs
        subIDs = find(model.S(:,i) < 0);
        % Calculating the carbon same as b4
        subCarbon = sum(abs(full(model.S(subIDs,i))).*metCarb(subIDs,1));
        bericht{i,1} ='IgG Heavy Chain Synthesis';
        constraint = abs(total_carbon_cons/subCarbon);
        if abs(minmax(i,2)) > constraint
            minmax_new(i,2) = constraint;
            mmAdj(i,2) = 1;
        end
        continue
    elseif (i == find(ismember(model.rxns,'igg_lc'))) %LC Synthesis
        % Get the coefficients for the biomass, not using the model_temp
        % because we want the coefpairs
        subIDs = find(model.S(:,i) < 0);
        % Calculating the carbon same as b4
        subCarbon = sum(abs(full(model.S(subIDs,i))).*metCarb(subIDs,1));
        bericht{i,1} ='IgG Light Chain Synthesis';
        constraint = abs(total_carbon_cons/subCarbon);
        if abs(minmax(i,2)) > constraint
            minmax_new(i,2) = constraint;
            mmAdj(i,2) = 1;  
        end
        continue
    end
    % evaluate if rxn is part of the measured constraints
    if ismember(i,varID_tmp)
        % If rxn is measured, then we don't set carbon constraints
        % as there are already lower and upper constraints to the reaction.
        bericht{i,1} = 'set_bnd'; 
        continue  
    % evaluate if rxn is part of the skip-rxns set       
    elseif ismember(i,skipRxns)
        bericht{i,1} = 'unconstr - oxPhos';
        continue
    % evaluate if rxn is balanced or unbalanced
    else
        % First, check if the reaction carries any flux under normal conditions
        if sum(abs(minmax(i,:))) == 0
            bericht{i,1} = 'no_flux';
            continue
        % Check wheather the reaction is aprt fo the drain set/exchange
        % between the cell and the enviroment.            
        elseif ismember(i,isDrain)    
            bericht{i,1} = 'balanced - drain';
            subIDs = find(model_tmp.S(:,i) < 0);
            prodIDs = find(model_tmp.S(:,i) > 0);

            % account for stoichiometric coefficients and add up carbon
            % As we set all the lower bounds for the exchange reactions to
            % 0, the atom constrain will take effect only on the upper
            % bound
            subCarbon = sum(abs(full(model_tmp.S(subIDs,i))).*metCarb(subIDs,1));
            prodCarbon = sum(abs(full(model_tmp.S(prodIDs,i))).*metCarb(prodIDs,1));

            if subCarbon == 0 && prodCarbon == 0 
                continue
            elseif subCarbon == 0
                constraint = abs(total_carbon_cons/prodCarbon);
            else
                constraint = abs(total_carbon_cons/subCarbon);
            end
            if model_tmp.rev(i,1) == 1
                if abs(minmax(i,1)) > constraint
                    minmax_new(i,1) = -constraint;
                    mmAdj(i,1) = 1;
                end
                if abs(minmax(i,2)) > constraint
                    minmax_new(i,2) = constraint;
                    mmAdj(i,2) = 1;
                end
            else
                if abs(minmax(i,2)) > constraint
                    if constraint > abs(minmax(i,1))
                        minmax_new(i,2) = constraint;
                        mmAdj(i,2) = 1;
                    end
                end
            end            
        else
            % Meaning, if the reaction is not 0 or drain
            
            % find substrates and products 
            subIDs = find(model_tmp.S(:,i) < 0);
            prodIDs = find(model_tmp.S(:,i) > 0);            
            % account for stoichiometric coefficients and add up carbon
            subCarbon = sum(abs(full(model_tmp.S(subIDs,i))).*metCarb(subIDs,1));
            prodCarbon = sum(abs(full(model_tmp.S(prodIDs,i))).*metCarb(prodIDs,1));
            
            if ismember(i,coa.coaRxnIDs)
                if coa.missMatch
                    warning('Stoichiometry of CoA is not matching. Smaller number is subtracted from carbon balance. Check reaction reaction manually.')
                    if coa.subAtom(ismember(coa.coaRxnIDs,i)) > coa.subCarbon(ismember(coa.coaRxnIDs,i))
                        subCarbon = subCarbon - coa.subAtom(ismember(coa.coaRxnIDs,i));
                        prodCarbon = prodCarbon - coa.prdAtom(ismember(coa.coaRxnIDs,i));
                    else
                        subCarbon = subCarbon - coa.subAtom(ismember(coa.coaRxnIDs,i));
                        prodCarbon = prodCarbon - coa.prdAtom(ismember(coa.coaRxnIDs,i));
                    end
                else
                    subCarbon = subCarbon - coa.subAtom(ismember(coa.coaRxnIDs,i));
                    prodCarbon = prodCarbon - coa.prdAtom(ismember(coa.coaRxnIDs,i));
                end
            end     
            
            % Account for the large molecules cytochromes, hemes, and quinones
            % and remove the atoms taking part there (need to check why
            % though???????)                    
            if ismember(i,lrgMlc.lmRxnIDs)
                subCarbon = subCarbon - (lrgMlc.sub(ismember(lrgMlc.lmRxnIDs,i)));
                prodCarbon = prodCarbon - lrgMlc.prd(ismember(lrgMlc.lmRxnIDs,i));
            end         
            
            % Check whether the number of atoms in the substrates or
            % products are equal to 0; even if the minmax is different tha
            % 0 the molecules could not have any carbons or nitrogens etc.                 
            if subCarbon == 0 && prodCarbon == 0
                bericht{i,1} = 'no C in the reaction';
                continue
            end

            % If in fact, the reaction is balanced then do the following:
            if subCarbon == prodCarbon     
                bericht{i,1} = 'balanced'; 
                constraint = abs(total_carbon_cons/subCarbon);
                if model_tmp.rev(i,1) == 1
                    if abs(minmax(i,1)) > constraint
                        minmax_new(i,1) = -constraint;
                        mmAdj(i,1) = 1;
                    end
                    if abs(minmax(i,2)) > constraint
                       minmax_new(i,2) = constraint;
                       mmAdj(i,2) = 1;
                    end
                else
                    if abs(minmax(i,2)) > constraint
                        if constraint > abs(minmax(i,1))
                            minmax_new(i,2) = constraint;
                            mmAdj(i,2) = 1;
                        end
                    end
                end     
            else
                
                % If the reaction is not balanced                
                % First check for the cofactor pairs              
                if find(ismember(cofMetIDs,find(model.S(:,i))))
                    % if the sum of all carbon is even then it is balanced 
                    % (not bulletproof but if off by only one carbon it works)
                    % atoms for the substrates, products and the coenzymes
                    % - cofactors
                    if mod(sum([subCarbon;prodCarbon;metCarb(cofMetIDs(ismember(cofMetIDs,find(model.S(:,i)))))]),2) == 0 
                        bericht{i,1} = 'balanced - CofPair';
                        % if balanced then its go with the smaller one (constraining less).
                        % This is because there seems to be only 1 cofactor in this equation 
                        % (e.g. atp -> amp, this is not picked up as a cofactor pair but atp 
                        % is a cofactor) so the carbon of atp/amp should not be considered when constraining
                        if subCarbon > prodCarbon 
                            constraint = abs(total_carbon_cons/prodCarbon);
                        else
                            constraint = abs(total_carbon_cons/subCarbon);
                        end
                        if model_tmp.rev(i,1) == 1
                            if abs(minmax(i,1)) > constraint
                                minmax_new(i,1) = -constraint;
                                mmAdj(i,1) = 1;
                            end
                            if abs(minmax(i,2)) > constraint
                               minmax_new(i,2) = constraint;
                               mmAdj(i,2) = 1;
                            end
                        else
                            if abs(minmax(i,2)) > constraint
                                if constraint > abs(minmax(i,1))
                                    minmax_new(i,2) = constraint;
                                    mmAdj(i,2) = 1;
                                end
                            end
                        end
                    else
                        
                        bericht{i,1} = 'unbalanced - CofPair';
                        % if unbalanced then go with the smaller one (constraining less)
                        % usually contains AMP which isn't considered to be
                        % a cofactor with no CoA though                        
                        if subCarbon < prodCarbon 
                            constraint = abs(total_carbon_cons/subCarbon);
                        else
                            constraint = abs(total_carbon_cons/prodCarbon);
                        end
                        if model_tmp.rev(i,1) == 1
                            if abs(minmax(i,1)) > constraint
                                minmax_new(i,1) = -constraint;
                                mmAdj(i,1) = 1;
                            end
                            if abs(minmax(i,2)) > constraint
                               minmax_new(i,2) = constraint;
                               mmAdj(i,2) = 1;
                            end
                        else
                            if abs(minmax(i,2)) > constraint
                                if constraint > abs(minmax(i,1))
                                    minmax_new(i,2) = constraint;
                                    mmAdj(i,2) = 1;
                                end
                            end
                        end
                        
                        % check whether there is a metabolite with X formula type.
                        flag = 0;    
                        f=1;
                        while (f<= size(subIDs,1) && flag == 0)
                            if (model.metFormulas{subIDs(f,1)} == 'X')
                                flag = 1;
                            end
                            f = f + 1;
                        end
                        f = 1;
                        while (f <= size(prodIDs,1) && flag == 0)
                            if (model.metFormulas{prodIDs(f,1)} == 'X')
                                flag = 1;
                            end
                            f = f + 1;
                        end
                        if (flag == 1)
                            bericht{i,1} = 'Contains X';
                            % Set as the constrain the larger one if there is an X
                            % metabolite; since normally the reaction would have
                            % been balanced - apprantly every time that's the case.
                            if subCarbon > prodCarbon
                                constraint = abs(total_carbon_cons/subCarbon);
                            else
                                constraint = abs(total_carbon_cons/prodCarbon);
                            end
                            if model_tmp.rev(i,1) == 1
                                if abs(minmax(i,1)) > constraint
                                    minmax_new(i,1) = -constraint;
                                    mmAdj(i,1) = 1;
                                end
                                if abs(minmax(i,2)) > constraint
                                    minmax_new(i,2) = constraint;
                                    mmAdj(i,2) = 1;
                                end
                            else
                                if abs(minmax(i,2)) > constraint
                                    if constraint > abs(minmax(i,1))
                                        minmax_new(i,2) = constraint;
                                        mmAdj(i,2) = 1;
                                    end
                                end
                            end
                        end
                    end
                else
                    % When the number of atoms in the substrate and product
                    % are differe, usually in rxns with metabolite X; and
                    % there is no coa or anything

                    % check whether there is a metabolite with X formula type.
                    flag = 0;    
                    f=1;
                    while (f<= size(subIDs,1) && flag == 0)
                        if (model.metFormulas{subIDs(f,1)} == 'X')
                            flag = 1;
                        end
                        f = f + 1;
                    end
                    f = 1;
                    while (f <= size(prodIDs,1) && flag == 0)
                        if (model.metFormulas{prodIDs(f,1)} == 'X')
                            flag = 1;
                        end
                        f = f + 1;
                    end
                    if (flag == 1)
                        bericht{i,1} = 'Contains X';
                        % Set as the constrain the larger one if there is an X
                        % metabolite; since normally the reaction would have
                        % been balanced - apprantly every time that's the case.
                        if subCarbon > prodCarbon
                            constraint = abs(total_carbon_cons/subCarbon);
                        else
                            constraint = abs(total_carbon_cons/prodCarbon);
                        end
                        if model_tmp.rev(i,1) == 1
                            if abs(minmax(i,1)) > constraint
                                minmax_new(i,1) = -constraint;
                                mmAdj(i,1) = 1;
                            end
                            if abs(minmax(i,2)) > constraint
                                minmax_new(i,2) = constraint;
                                mmAdj(i,2) = 1;
                            end
                        else
                            if abs(minmax(i,2)) > constraint
                                if constraint > abs(minmax(i,1))
                                    minmax_new(i,2) = constraint;
                                    mmAdj(i,2) = 1;
                                end
                            end
                        end
                    else
                         % When there is no CofPair or X (there is probably
                         % EPO or HX but we didn't account for it yet :(
                        bericht{i,1} = 'unbalanced';
                        if subCarbon > prodCarbon
                            constraint = abs(total_carbon_cons/subCarbon);
                        else
                            constraint = abs(total_carbon_cons/prodCarbon);
                        end
                        if model_tmp.rev(i,1) == 1
                            if abs(minmax(i,1)) > constraint
                                minmax_new(i,1) = -constraint;
                                mmAdj(i,1) = 1;
                            end
                            if abs(minmax(i,2)) > constraint
                                minmax_new(i,2) = constraint;
                                mmAdj(i,2) = 1;
                            end
                        else
                            if abs(minmax(i,2)) > constraint
                                if constraint > abs(minmax(i,1))
                                    minmax_new(i,2) = constraint;
                                    mmAdj(i,2) = 1;
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    


end

%% qunatify how much solution space was constraint compared to original minmax
mm_diff = minmax(:,2)-minmax(:,1);
mm_new_diff = minmax_new(:,2)-minmax_new(:,1);
diff = 100*(1-(sum(mm_new_diff)/sum(mm_diff)));
        
%% variable complier
unbRxns = sum(ismember(bericht,'unbalanced'))+sum(ismember(bericht,'unbalanced - CofPair'));
balRxns = sum(ismember(bericht,'balanced'))+sum(ismember(bericht,'balanced - X'))...
    +sum(ismember(bericht,'balanced - drain'))+sum(ismember(bericht,'unconstr - oxPhos'));
setRxns = sum(ismember(bericht,'set_bnd'));
nfRxns = sum(ismember(bericht,'no_flux'));
Unknowns = sum(ismember(bericht,'Contains X'));

carbonConst.title = '*** carbon closure constraining ***';
carbonConst.dateAndTime = datetime();
carbonConst.constraintsName = condition;
if isfield(model,'description')
    carbonConst.model = model.description;
end
carbonConst.totalCarbonConsumed = total_carbon_cons;
if exist('actual_total_carbon_cons','var')
    if ~isempty(actual_total_carbon_cons)
        carbonConst.actualTotalCarbonConsumed = actual_total_carbon_cons;
    end
end
carbonConst.report.summary = bericht;
carbonConst.report.unbalanced = unbRxns;
carbonConst.report.balanced = balRxns;
carbonConst.report.setRxns = setRxns;
carbonConst.report.noFlux = nfRxns;
carbonConst.report.Unknowns = Unknowns;
carbonConst.bndAjusted.rawData = mmAdj;
carbonConst.bndAjusted.numberOf_LB_adj = sum(mmAdj(:,1));
carbonConst.bndAjusted.numberOf_UB_adj = sum(mmAdj(:,2));
carbonConst.CoAevaluation = coa;
carbonConst.carbonCount = carbonCount;
carbonConst.cofactorPairs = cofPairs;
carbonConst.originalMinMax = minmax;
carbonConst.newMinMax = minmax_new;
carbonConst.solSpaceReduction = strcat(num2str(diff),'%');

%% Test if model solves under new conditions

model.lb = minmax_new(:,1);
model.ub = minmax_new(:,2);

% for i = 1 : 3561
%     model.lb(1:i) = minmax_new(1:i,1);
%     model.ub(1:i) = minmax_new(1:i,2);
%     sol = optimizeCbModel(model,'max')
% if isnumeric(sol.origStat)
%     if sol.origStat < 1
%         error('Model does not solve. Think about relaxing the bounds')
%     else
%         disp(strcat('model solves. sol.x =',' ',num2str(sol.f)))
%     end
% else
%     if ismember(cellstr(sol.origStat),{'INFEASIBLE'})
%         warning('Model does not solve. Think about relaxing the bounds')
%     else
%         disp(strcat('model solves. sol.x =',' ',num2str(sol.f)))
%     end
% end
% end

sol = optimizeCbModel(model,'max');

if isnumeric(sol.origStat)
    if sol.origStat < 1
        error('Model does not solve. Think about relaxing the bounds')
    else
        disp(strcat('model solves. sol.x =',' ',num2str(sol.f)))
    end
else
    if ismember(cellstr(sol.origStat),{'INFEASIBLE'})
        warning('Model does not solve. Think about relaxing the bounds')
    else
        disp(strcat('model solves. sol.x =',' ',num2str(sol.f)))
    end
end

carbonConst.solutionNew = sol;

%% constrain uptakes based on passed in limit (e.g. lost/closed carbon)

if exist('limEX','var')
    if ~isempty(limEX)
        carbonConst.EXrxnAdjust = '*** EX rxns adjustment ***';
        carbonConst.limEX = limEX;
        
        % get all the max values for the EX rxns not set
        varMaxN = minmax_new(isDrainNew,2);
        % get met indices
        metID_ind = find(model.S(:,isDrainNew));
        % get metIDs from Met indices
        [metID_EX,~] = ind2sub(size(model.S),metID_ind);
        % get rxnIDs where rxn actually has a flux and ignore any fluxes
        % that are set to 0
        tmpRxnID = isDrainNew(find(varMaxN));
        % how many rxns have been set to 0?
        numEXzero = size(isDrainNew,1)-size(tmpRxnID,1);
        % get metIDs from EX rxn IDs so carbon per met can be determined
        tmpMetID = metID_EX(find(varMaxN));
        % calculate new bnds
        limEXnew = limEX./metCarb(tmpMetID);
        % get rid of Inf
        limEXnew(~isfinite(limEXnew)) = 0;
        % find which reactions needf adjustment
        adjEX = tmpRxnID(limEXnew < varMaxN(ismember(isDrainNew,tmpRxnID)));
        %  check if there is any bnds to be adjusted
        if ~isempty(adjEX)
            % make calculations as to how much solution space shrinks
            % compared to original
            diffLim = 100*(1-(sum(limEXnew(ismember(tmpRxnID,adjEX)))/sum(minmax_new(adjEX,2))));
            limEXadj = limEXnew(ismember(tmpRxnID,adjEX));
            minmax_new(adjEX,2) = limEXnew(ismember(tmpRxnID,adjEX));
            
            mm_lim_diff = minmax_new(:,2)-minmax_new(:,1);
            diff_n = 100*(1-(sum(mm_lim_diff)/sum(mm_diff)));
            
            model.lb = minmax_new(:,1);
            model.ub = minmax_new(:,2);
            sol = optimizeCbModel(model,'max');

            if isempty(sol.f)
                warning('model does not solve. think about relaxing the bounds')
            else
                disp(strcat('model solves. sol.x =',' ',num2str(sol.f)))
            end
            
            carbonConst.numAdjusted = size(adjEX,1);
            carbonConst.solSpaceReductionEX  = strcat(num2str(diff_n),'%');
            carbonConst.EXspaceReduction = strcat(num2str(diffLim),'%');
        else
            warning('no adjustment of EXrxns was done as already constrained more than limEX')
            carbonConst.warning = 'no adjustment was done';
        end
    end
end

%% run minmax on new minmax bnds and see if solution space is reduced further
if runFinalMinMax
    model.lb = minmax_new(:,1);
    model.ub = minmax_new(:,2);
    minmax_post_new = runMinMax_GF(model);
    carbonConst.postNmm = '*** post new minmax ***';
    carbonConst.postNewMinMax = minmax_post_new;
    % qunatify how much solution space was constraint compared to original minmax
    mm_post_diff = minmax_post_new(:,2)-minmax_post_new(:,1);
    diff = 100*(1-(sum(mm_post_diff)/sum(mm_diff)));
    carbonConst.solSpaceReductionPostCC = strcat(num2str(diff),'%');

    % test if new minmax solves
    model.lb = minmax_post_new(:,1);
    model.ub = minmax_post_new(:,2);
    sol = optimizeCbModel(model,'max');

    if isempty(sol.f)
        warning('model does not solve. think about relaxing the bounds')
    else
        disp(strcat('model solves. sol.x =',' ',num2str(sol.f)))
    end

    carbonConst.solutionPostNew = sol;
end

