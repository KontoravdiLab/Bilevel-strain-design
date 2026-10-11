% The lower and upper bounds for each of the reactions need to be set to a
% value larger than the total carbon/nitrogen or smth uptake, as only then
% the bounds will bo constrained to a realistic value.

%% Parameter set up and clear
clear; clc; close all;
changeCobraSolver('gurobi')
max_bound = 100;

dynamic = 0;
doctor(1) = {'SUE'}; % Susie FeedC
doctor(1) = {'SAR_FeedC'}; % Sarantos FeedC
doctor(3) = {'SAR_FeedAll'}; % Sarantos FeedAll
doctor(4) = {'SAR_FeedAll40'}; % Sarantos FeedAll plus 40%
% doctor(5) = {'WID'}; % Combination of the previous set of data

% dynamic = 1;
% doctor(1) = {'SUE'};
% doctor(2) = {'SAR_FeedC'};
% doctor(3) = {'SAR_FeedAll'};
% doctor(4) = {'SAR_FeedAll40'};

% doctor = {'ROD'};

%% Model set-up

% MODEL = 'GENOME_SCALE';
% load Cho_GeM_v4.mat

% MODEL = 'GENOME_SCALE_NoEPO';
% load Cho_GeM_NoEpo.mat;

% MODEL = 'K1';
% load Cho_GeM_K1.mat; model = modelnew;
% model.rev = checkReversability(model);

% MODEL = 'K1_NoEPO';
% load Cho_GeM_K1_NoEpo.mat;

MODEL = 'Secretome';
load iCHO2291_secretome_curated.mat;
model.rev = checkReversability(model);


% MODEL = 'DG_44';
% load Cho_GeM_DG44.mat; model = modelnew;model.rev = checkReversability(model);
% load Cho_DG44.mat

% MODEL = 'S';
% load Cho_GeM_S.mat; model = modelnew;

% MODEL = 'Core';
% % load e_coli_core.mat; model = e_coli_core; model.csense = model.csense';

% MODEL = 'iJO1366';
% load iJO1366.mat; model = iJO1366; model.csense = model.csense';

%% Atom Count declaration
atomCount = findMetCarbon(model);

wmax = size(doctor,2);
for ww = 1:wmax
    doc = doctor{ww};

    if dynamic == 1
        
        %%%%% Reset the model bounds
        rev = checkReversability(model);
        model.rev = rev;
        [model,~] = resetModel(model,-max_bound,max_bound,model.rev);
        model.lb(model.ExchRxnBool) = 0;
        model.lb(model.DMRxnBool) = 0;
        model.lb(model.SinkRxnBool) = -max_bound;
        model.ub(model.ExchRxnBool) = max_bound;
        model.ub(model.DMRxnBool) = max_bound;
        model.ub(model.SinkRxnBool) = max_bound;
        
        [rxn,bounds] = Sarantos_Bounds_14Day(doc);
%         [rxn,bounds] = Sarantos_Bounds_14Day('All40');

        days = size(bounds,2);
        for stepNo = 1 : days - 1
            bnds_rxns(:,1) = rxn;
            lb = min(bounds(:,stepNo:stepNo+1),[],(2));
            ub = max(bounds(:,stepNo:stepNo+1),[],(2)); ub(24:35) = 10;                    
            for j=1:size(bounds,1)                
                if ub(j) <= lb(j)
                    ub(j) = max_bound;
                end
                bnds_rxns(j,2) = {lb(j)};
                bnds_rxns(j,3) = {ub(j)};
            end
            model_temp = set_bounds(model,lb,ub,rxn);

%             model_temp.c = 0 * model_temp.c;
%             model_temp.c(4005) = 1;
            sol = optimizeCbModel(model_temp);
            sol.f
            if sol.stat == 0 || sol.obj == 0
                continue
%             else
%                 continue
            end
            condition  = 'Dynamic';        
            smallMets = findSmallMets(model);
            metPairs = countMetPairs_max(model,smallMets.smallMetIDs,true,'[');
            cofPairs = findCofactorPairs_max(metPairs);
            minmax = runMinMax_GF(model_temp); 
            minmax = fixMinMax(minmax);  
            
            runFinalMinMax = false;        
            skipRxns = [];
            
            [carbonConst] = carbonConstraints(model_temp,bnds_rxns,condition,atomCount,cofPairs,minmax,runFinalMinMax,[],[],skipRxns);
            [nitrogenConst] = nitrogenConstraints(model_temp,bnds_rxns,condition,atomCount,cofPairs,minmax,runFinalMinMax,[],[],skipRxns);
            
            
            New_MinMax = selectMinMax(carbonConst, nitrogenConst);
            mm = New_MinMax;
            
            for i = 1:size(mm)
                A_MM = mm{i,1};
                B_MM = mm{i,2};
                if (iscell(A_MM))
                    A_MM = A_MM{1,1};
                end
                if (iscell(B_MM))
                    B_MM = B_MM{1,1};
                end
                Bounds(i,1) = A_MM;
                Bounds(i,2) = B_MM;
            end
            
            model_1 = model_temp;
            model_1.lb = Bounds(:,1);
            model_1.ub = Bounds(:,2);
            
            sol_temp = optimizeCbModel(model_1);
            if sol_temp.f == 0
                continue
            end
            FVA = runMinMax_GF(model_1);            
            
            fva = sprintf('FVA_EX_%s_%s_DAY_%d.mat',doc,MODEL,stepNo);
%             save(fva,'FVA')  
        
            %% Save the newMinMax in the folder
            % > > > > > > > > > SAVING  WORKSPACE > > > > > > > > > > > > > > > > > > > > > > > > > > > > >
%             [dateStr, timeStr] = getDateTimeStrings(date,clock);     
%             eval(['save ./','/',dateStr,'_',timeStr,'_',mfilename,'_',doc,'_',MODEL,'_DAY_',stepNo,'.mat;'])
%             fname = sprintf('%s_%s_%s_%s_%s_DAY_%d.mat',dateStr,timeStr,mfilename,doc,MODEL,stepNo);
%             save(fname)
            % < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < <
        
        end
        
    elseif dynamic == 0 
        
        [lbg,ubg,rxng,modelg] = FBA_Bounds(model,1,max_bound,doc);
        lbg = lbg';
        ubg = ubg';
        rxng = rxng';
        bndsg = [rxng lbg ubg];
        modelg = set_bounds(modelg,lbg,ubg,rxng);

        
        for phase = 1:1
%             d_phase = sprintf('%d',phase);
%             doc = strcat(doc,d_phase);
%             [lbg,ubg,rxng,modelg] = RODS_Boudnds(model,phase,max_bound);
%             lbg = lbg;
%             ubg = ubg;
%             rxng = rxng';
%             bndsg = [rxng lbg ubg];
%             modelg = set_bounds(modelg,lbg,ubg,rxng);

%         sol1 = optimizeCbModel(modelg);
%         modelg.ub(find(ismember(modelg.rxns,'EX_fru_e_'))) = -0.001; % EX_cys_L_e_
%         modelg.lb(find(ismember(modelg.rxns,'EX_fru_e_'))) = -0.001; % EX_cys_L_e_
%         sol2 = optimizeCbModel(modelg);
%         modelg.ub(find(ismember(modelg.rxns,'EX_fru_e_'))) = -0.01; % EX_cys_L_e_
%         modelg.lb(find(ismember(modelg.rxns,'EX_fru_e_'))) = -0.05; % EX_cys_L_e_
%         sol3 = optimizeCbModel(modelg);
%         modelg.ub(find(ismember(modelg.rxns,'EX_fru_e_'))) = -0.01; % EX_cys_L_e_
%         modelg.lb(find(ismember(modelg.rxns,'EX_fru_e_'))) = -0.1; % EX_cys_L_e_
%         modelg.lb(find(ismember(modelg.rxns,'G6PPer'))) = 0.01; % EX_cys_L_e_
% 
%         sol4 = optimizeCbModel(modelg);
%         CCC(modelg,({'fru[e]';'fru[c]';'g6p[c]';'g6p[r]';'glc_D[r]';'udpg[c]';'gdpmann[c]';'udpg[c]';}),sol1.x,sol2.x,sol3.x,sol4.x)

        
        % Additional Constraints
%         modelg = changeRxnBounds(modelg, 'G6PDH1rer', 0, 'l'); % PPP 1 irriversible cycles; THIS ONE
        % modelg = changeRxnBounds(modelg, 'G6PDH2rer', 0, 'l'); % PPP 1 irriversible cycles; THIS ONE
        % modelg = changeRxnBounds(modelg, 'G6PDH2r', 0, 'l'); % PPP 1 irriversible cycles; THIS ONE
        % modelg = changeRxnBounds(modelg, 'GND', 0, 'l'); % PPP 3 irriversible cycles; THIS ONE
        % modelg = changeRxnBounds(modelg, 'GNDer', 0, 'l'); % PPP 3 irriversible cycles; THIS ONE
        % modelg = changeRxnBounds(modelg, 'ATPS4m_cho', -max, 'l'); % ATP Synthase is reversible; THIS ONE
        
        %%% Stationary Phase 
%         [lbs,ubs,rxns,models] = FBA_Bounds(model,2,max_bound,doc);
%         lbs = lbs';
%         ubs = ubs';
%         rxns = rxns';
%         bndss = [rxns lbs ubs];
%         models = set_bounds(models,lbs,ubs,rxns);
        
        % Additional Constraints
        % models = changeRxnBounds(models, 'G6PDH1rer', 0, 'l'); % PPP 1 irriversible cycles; THIS ONE
        % models = changeRxnBounds(models, 'G6PDH2rer', 0, 'l'); % PPP 1 irriversible cycles; THIS ONE
        % models = changeRxnBounds(models, 'G6PDH2r', 0, 'l'); % PPP 1 irriversible cycles; THIS ONE
        % models = changeRxnBounds(models, 'GND', 0, 'l'); % PPP 3 irriversible cycles; THIS ONE
        % models = changeRxnBounds(models, 'GNDer', 0, 'l'); % PPP 3 irriversible cycles; THIS ONE
        % models = changeRxnBounds(models, 'ATPS4m_cho', -max, 'l'); % ATP Synthase is reversible; THIS ONE
        
        %% Condition declaration
        condition1  = 'Growth Phase' ;
%         condition2  = 'Stationary Phase' ;
        
        %% Cofactor pairs declaration    
        smallMets = findSmallMets(model);
        % IMPORTANT, for our cho compartmets use [ not _ !!!
        metPairs = countMetPairs_max(model,smallMets.smallMetIDs,true,'[');
        cofPairs = findCofactorPairs_max(metPairs);
        
        %% MinMax declaration
        minmaxg = runMinMax_GF(modelg); 
        minmaxg = fixMinMax(minmaxg);  
%         minmaxs = runMinMax_GF(models); 
%         minmaxs = fixMinMax(minmaxs);  
   
        %% Run final MinMax declaration
        % False because we run a minmax in the end
        runFinalMinMax = false;
        
        %% Declare which reactions to skip
        skipRxnsIDs = { 'CYOOm2_cho'...
        'CYOOm3_cho'...
        'CYOR_u10m_cho'...
        'NADH2_u10m_cho'...
        'PPA'...
        'PPAm'...
        'ATPS4m_cho'...
        'FADH2ETC'...
        'GLYC3PFADm'...
         }';
        
        c = 1;
        
        for i = 1:size(skipRxnsIDs,1)
            for j = 1:size(model.rxns,1)
                if strcmp(skipRxnsIDs(i),model.rxns(j))
                    skipRxns(c) = j;
                    c = c + 1;
                end
            end
        end
%         skipRxns = [];    
        
        %% Carbon Constrain Algorithm
        [carbonConstExp] = carbonConstraints(modelg,bndsg,condition1,atomCount,cofPairs,minmaxg,runFinalMinMax,[],[],skipRxns);
%         [carbonConstStat] = carbonConstraints(models,bndss,condition2,atomCount,cofPairs,minmaxs,runFinalMinMax,[],[],skipRxns);
        
        %% Nitrogen Constrain Algorithm 
        [nitrogenConstExp] = nitrogenConstraints(modelg,bndsg,condition1,atomCount,cofPairs,minmaxg,runFinalMinMax,[],[],skipRxns);
%         [nitrogenConstStat] = nitrogenConstraints(models,bndss,condition2,atomCount,cofPairs,minmaxs,runFinalMinMax,[],[],skipRxns);
        
        %% Combine the two minmaxes
        % Experimental MinMax
        New_Exp_MinMax = selectMinMax(carbonConstExp, nitrogenConstExp);
        % Stationary MinMax
%         New_Stat_MinMax = selectMinMax(carbonConstStat, nitrogenConstStat);
        
        %% Measure the combined solution space reduction
%         for i = 1: size(New_Exp_MinMax,1)
%             A_Exp = New_Exp_MinMax{i,2};
%             B_Stat = New_Exp_MinMax{i,1};
%             if (iscell(B_Stat))
%                 B_Stat = B_Stat{1,1};
%             end
%             if (iscell(A_Exp))
%                 A_Exp = A_Exp{1,1};
%             end
%             mm_new_diff(i) = A_Exp - B_Stat;
%             A_Exp = New_Stat_MinMax{i,2};
%             B_Stat = New_Stat_MinMax{i,1};
%             if (iscell(B_Stat))
%                 B_Stat = B_Stat{1,1};
%             end
%             if (iscell(A_Exp))
%                 A_Exp = A_Exp{1,1};
%             end
%             mm_new_diff(i) = A_Exp - B_Stat;
%         end
%         mm_diff = minmaxg(:,2)-minmaxg(:,1);
%         diff_exp = 100*(1-(sum(mm_new_diff)/sum(mm_diff)))
        
%         mm_diff = minmaxs(:,2)-minmaxs(:,1);    
%         diff_stat = 100*(1-(sum(mm_new_diff)/sum(mm_diff)))
        
        %% Convert the minmax from Cell to Double
        mm_exp = New_Exp_MinMax;
%         mm_stat = New_Stat_MinMax;
        
        for i = 1:size(mm_exp)
            A_Exp = mm_exp{i,1};
            B_Exp = mm_exp{i,2};
%             A_Stat = mm_stat{i,1};
%             B_Stat = mm_stat{i,2};
            if (iscell(A_Exp))
                A_Exp = A_Exp{1,1};
            end
            if (iscell(B_Exp))
                B_Exp = B_Exp{1,1};
            end
%             if (iscell(A_Stat))
%                 A_Stat = A_Stat{1,1};
%             end
%             if (iscell(B_Stat))
%                 B_Stat = B_Stat{1,1};
%             end
            Exponential_Bounds(i,1) = A_Exp;
            Exponential_Bounds(i,2) = B_Exp;
%             Stationary_Bounds(i,1) = A_Stat;
%             Stationary_Bounds(i,2) = B_Stat;
        end
        
        %% Run a new minmax to constrain the combined bounds - carbon and nitrogen
        model_1 = modelg;
        model_1.lb = Exponential_Bounds(:,1);
        model_1.ub = Exponential_Bounds(:,2);
        
        sol = optimizeCbModel(model_1);

        
%         model_2 = models;
%         model_2.lb = Stationary_Bounds(:,1);
%         model_2.ub = Stationary_Bounds(:,2);
        
        FVA_EX = runMinMax_GF(model_1);
        FVA_EX = fixMinMax(FVA_EX);
%         FVA_ST = runMinMax_GF(model_2);


%% FVA TEST

        model_1.lb = FVA_EX(:,1);
        model_1.ub = FVA_EX(:,2);
        sol = optimizeCbModel(model_1);
        if sol.stat ~= 1
            1
        end
%     for i =2480:3252
%         model_1.lb(i) = FVA_EX(i,1);
%         model_1.ub(i) = FVA_EX(i,2);
%         sol = optimizeCbModel(model_1);
%         if sol.stat ~= 1
%             1
%         end
%     end
        
        %% Calculate the total difference - combined vs original
%         for i = 1: size(FVA_EX,1)
%             A_1 = FVA_EX(i,2);
%             B_1 = FVA_EX(i,1);
%             mm_new_diff_exp(i) = A_1 - B_1;
%             A_1 = FVA_ST(i,2);
%             B_1 = FVA_ST(i,1);
%             mm_new_diff_stat(i) = A_1 - B_1;
%         end
%         mm_diff = minmaxg(:,2)-minmaxg(:,1);
%         diff_exp = 100*(1-(sum(mm_new_diff_exp)/sum(mm_diff)));
%         mm_diff = minmaxs(:,2)-minmaxs(:,1);
%         diff_stat = 100*(1-(sum(mm_new_diff_stat)/sum(mm_diff)));

        fva_ex = sprintf('FVA_EX_%s_%s_NEW.mat',doc,MODEL);
%         fva_st = sprintf('FVA_ST_%s_%s.mat',doc,MODEL);
        save(fva_ex,'FVA_EX')
%         save(fva_st,'FVA_ST')
        
        %% Save the newMinMax in the folder
        % > > > > > > > > > SAVING  WORKSPACE > > > > > > > > > > > > > > > > > > > > > > > > > > > > >
        [dateStr, timeStr] = getDateTimeStrings(date,clock);     
        eval(['save ./','/',dateStr,'_',timeStr,'_',mfilename,'_',doc,'_',MODEL,'.mat;'])
        % < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < <
        
        end
    end
    
end

%% Test Bounds
% a = modelg;
% for i =1:size(FVA_EX,1)
%     a.lb(i) = FVA_EX(i,1);
%     x = optimizeCbModel(a);
%     if(x.f <0.1)
%         Exponential_Bounds(i,1)
%         FVA_EX(i,1)
%         i
%         pause
%     end
% end



%% Plot the results of FVA-Carbon-etc
carbon_nitrogen_plot(minmaxg,carbonConstExp,nitrogenConstExp,model)