function [essential_rxns_ko, essential_rxns_up, essential_rxns_down] = GSA_Genes(model,genes_up,genes_ko,genes_down,up,down)

%% Initialization of Gene sets
genes_up_init = unique(genes_up);
genes_ko_init = unique(genes_ko);
genes_down_init = unique(genes_down);

if isnumeric(genes_ko(1))
    genes = [unique(genes_up_init); unique(genes_ko_init); unique(genes_down_init)];
    num_up = size(unique(genes_up_init),1);
    num_null = size(unique(genes_ko_init),1);
else
    genes = [unique(genes_up_init) unique(genes_ko_init)];
    num_up = size(unique(genes_up_init),2);
end
imax = size(genes,1);

%% Get the initial FBA solution with all the regulations
DM_IgG = findRxnIDs(model,'DM_igg[g]');
bounds_init_wt(:,1) = model.lb;
bounds_init_wt(:,2) = model.ub;
bounds_init_es = bounds_init_wt;
bounds_init_es(genes_up_init,1) = bounds_init_wt(genes_up_init,1) + (bounds_init_wt(genes_up_init,2)-bounds_init_wt(genes_up_init,1)) * up;
bounds_init_es(genes_ko_init,:) = 0;
bounds_init_es(genes_down_init,2) = bounds_init_wt(genes_down_init,2) - (bounds_init_wt(genes_down_init,2)-bounds_init_wt(genes_down_init,1)) * down;
model.lb = bounds_init_es(:,1);
model.ub = bounds_init_es(:,2);
sol_init = optimizeCbModel(model);
bio_init = sol_init.f;
dm_init = sol_init.x(DM_IgG);

%% Detect the  genetic regulations that have no effect in the biomass or IgG production. 
counter = 0;
genes_up_init_es = [];
genes_ko_init_es = [];
genes_down_init_es = [];

for i = 1:imax
    genes_up_init_new = genes_up_init;
    genes_ko_init_new = genes_ko_init;
    genes_down_init_new = genes_down_init;

    is_up = false;
    is_down = false;
   
    if (i <= num_up)
        genes_up_init_new = genes_up_init(1:end ~= i);
        is_up = true;
    elseif (i <= num_up + num_null)
        genes_ko_init_new = genes_ko_init(1:end ~= (i - num_up));
    else
        genes_down_init_new = genes_down_init(1:end ~= (i - num_up - num_null));
        is_down = true;
    end
    
    % New optimisation
    bounds_init_es = bounds_init_wt;
    bounds_init_es(genes_up_init_new,1) = bounds_init_wt(genes_up_init_new,1) + (bounds_init_wt(genes_up_init_new,2)-bounds_init_wt(genes_up_init_new,1)) * up;
    bounds_init_es(genes_ko_init_new,:) = 0;    
    bounds_init_es(genes_down_init_new,2) = bounds_init_wt(genes_down_init_new,2) - (bounds_init_wt(genes_down_init_new,2)-bounds_init_wt(genes_down_init_new,1)) * down;
    model.lb = bounds_init_es(:,1);
    model.ub = bounds_init_es(:,2);
    sol = optimizeCbModel(model);  
    sol.f;
    sol.x(DM_IgG);
    
    accurancy_level = 0.9999;
    if ( sol.x(DM_IgG)/dm_init > accurancy_level && sol.f/bio_init > accurancy_level ) 
        counter = counter + 1;
        if is_up 
            genes_up_init_es(counter,1) = genes(i);
        elseif is_down
            genes_down_init_es(counter,1) = genes(i);
        else
            if ~exist('genes_up_init_es','var') || ~exist('genes_down_init_es','var')
                genes_ko_init_es(counter,1) = genes(i);
            else
                genes_ko_init_es(counter - size(genes_up_init_es,1) - size(genes_down_init_es,1),1) = genes(i);
            end
        end 
    end    
end

%% Collect the essential reactions
if exist('genes_up_init_es','var')
    essential_rxns_up = genes_up_init(~ismember(genes_up_init,genes_up_init_es));
else
    essential_rxns_up = genes_up_init;
end
if exist('genes_ko_init_es','var')
    essential_rxns_ko = genes_ko_init(~ismember(genes_ko_init,genes_ko_init_es));
else
    essential_rxns_ko = genes_ko_init;
end
if exist('genes_down_init_es','var')
    essential_rxns_down = genes_down_init(~ismember(genes_down_init,genes_down_init_es));
else
    essential_rxns_down = genes_down_init;
end

end