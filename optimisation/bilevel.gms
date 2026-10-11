*****************************************************************************************
*  bilevel.gms  -  bilevel strain design for IgG production in CHO-K1
*
*  REPRESENTATIVE SCRIPT. The published scenarios come from many runs of this
*  program, each run with its own settings. Before running, set the options in
*  the RUN SETTINGS block below: dataset, algorithm, the e_up and k ranges and
*  the time limit. The ranges used across the runs are given next to each one.
*  Every combination listed there was run (over 1,300 optimisation runs).
*
*  Model    iCHO_K1_used (published CHO-K1 reconstruction, 700 blocked and 71
*           EPO reactions removed; 3952 reactions), split into irreversible
*           reactions. Biomass = R3930 (biomass_cho_producing),
*           IgG = E3934 (DM_igg[g]).
*  Inputs   inputs/<dataset folder>/  (stoichiometry, bounds, reaction types)
*  Run      from this folder:   gams bilevel.gms
*  Solver   CPLEX 12.9.0.0 under GAMS 28.2.0 (commercial licence required)
*  Output   RXNOptimization_<dataset>_<algorithm>.txt, one block per (e_up, k)
*           with up to five alternative regulation sets, written as MATLAB
*           index assignments (r_up, r_null) for evaluation by FBA.
*
*  Formulation: the inner FBA problem is replaced by its dual and joined to
*  the outer problem by a strong-duality constraint (Primal_dual); products of
*  dual variables and binaries are linearised with big-M constraints (M = 100).
*
*  Reversibility constraint (ForBackRxn): this file uses the corrected form,
*  which links the forward and backward halves of each reversible reaction.
*  Earlier runs used a form that did not link them. Every reported scenario was
*  re-evaluated by FBA, so its reported values do not depend on this constraint.
*
*  This file follows the scripts run for D1, D3 and D4. The D2 runs differed:
*  Bio_Igg also excluded reaction R755 from regulation (added to constraint_flux)
*  and used CPLEX option file 4; Igg_Igg applied dual1b/dual1c to the IgG
*  reaction (E3934) instead of the biomass reaction.
*****************************************************************************************

*================================== RUN SETTINGS =======================================

* Dataset: D1 = SUE (CHO-T, FeedC), D2 = SAR_FeedC (CHO-GS46, FeedC),
*          D3 = SAR_FeedAll (CHO-GS46, FeedAll), D4 = SAR_FeedAll40 (FeedAll + 40%)
$set dataset   D4

* Algorithm: Bio_Igg = inner problem maximises biomass (growth-coupled design)
*            Igg_Igg = inner problem maximises IgG
$set algorithm Bio_Igg

* Upregulation coefficient e_up: an upregulated reaction must carry at least
* LB + e_up*(UB - LB). Runs used 0.4 to 1.1 in steps of 0.1
* (a few Igg_Igg runs also at 0.65 and 0.95).
$set eup_from  0.4
$set eup_to    1.1
$set eup_step  0.1

* Number of regulations k (knockouts + upregulations). Runs used 10 to 225:
* 10 to 150 in steps of 10, and 25 to 225 in steps of 25.
$set k_from    25
$set k_to      225
$set k_step    25

* Downregulation coefficient e_down. Downregulation is disabled (y_low fixed to
* 0 below), so e_down has no effect. It was 0.3 in all runs; one exploratory
* Igg_Igg batch enabled downregulation with e_down = 0.3, 0.6 and 0.9.
$set edown     0.3

* Solver time limit per solve, seconds. Runs used 10000 to 100000.
$set reslim    10000

*=======================================================================================

$ifThen %dataset% == D1
$set myroot    inputs/Susie/
$set doc       SUE
$set maxgrowth 0.012345
$elseIf %dataset% == D2
$set myroot    inputs/SarantosFeedC/
$set doc       SAR_FeedC
$set maxgrowth 0.029705
$elseIf %dataset% == D3
$set myroot    inputs/SarantosFeedAll/
$set doc       SAR_FeedAll
$set maxgrowth 0.029705
$elseIf %dataset% == D4
$set myroot    inputs/SarantosFeedAll40/
$set doc       SAR_FeedAll40
$set maxgrowth 0.029705
$else
$abort dataset must be D1, D2, D3 or D4
$endIf

$ifThen not %algorithm% == Bio_Igg
$ifThen not %algorithm% == Igg_Igg
$abort algorithm must be Bio_Igg or Igg_Igg
$endIf
$endIf

Options
* minimise the size of the LST file
        limrow = 0
        limcol = 0
        solprint = off
        sysout = off
* solution and solver options
        iterlim = 1000000
        decimals = 8
        reslim = %reslim%
        work = 50000000;

******************************** SET DEFINITIONS ****************************************
*       i                  = metabolites
*       j                  = reactions (irreversible split)
*       exchange(j)        = uptake reactions in the medium
*       Biom(j)            = biomass reaction
*       constraint_flux(j) = reactions excluded from regulation
*       NoGenes(j)         = reactions without a gene association
*       index              = store for alternative solutions
*****************************************************************************************
Sets

i
$include "%myroot%ModelMETS.txt"

j
$include "%myroot%ModelRXNames.txt"

exchange(j)
$include "%myroot%ModelUptakes.txt"

Biom(j) / 'R3930' /

constraint_flux(j)
/
* Biomass
'R3930'
$ifThen %algorithm% == Bio_Igg
* IgG production
'E3934'
$endIf
/

NoGenes(j)
$include "%myroot%ModelNoGenes.txt"

index /1*100/
;

Alias(j,jj);

****************************** PARAMETERS ***********************************************
*       s(i,j)      = stoichiometric coefficient of metabolite i in reaction j
*       rxntype(j)  = 0 irreversible, 1 forward, 2 backward, 4 exchange
*       Reverse(j)  = index of the partner half of a reversible reaction
*       LB(j)/UB(j) = flux bounds for the dataset
*****************************************************************************************
Parameters

s(i,j)
$include "%myroot%ModelS.txt"

rxntype(j)
$include "%myroot%ModelRXNStype.txt"

Reverse(j)
$include "%myroot%ModelBackwardRxns.txt"

LB(j) Lower bounds containing the medium and the feed
$include "%myroot%ModelLowerBound_%doc%.txt"

UB(j) Upper bounds containing the medium and the feed
$include "%myroot%ModelUpperBound_%doc%.txt"

maxb, maxi, cb, ci
e_up, e_down

* interventions already found, for the integer cuts
done_up(index,j)
done_ko(index,j)
done_down(index,j)

* counters used when writing the output
w1,w2,w3
jlnj

null_par(j), up_par(j), down_par(j)
;

***************************** VARIABLES *************************************************
*       v(j)                 = flux
*       z                    = outer objective (IgG flux)
*       lambda(i)            = duals of the stoichiometric constraints
*       muL(j)/muU(j)        = duals of the flux bounds
*       theta(j)/phi(j)      = duals of the up- and downregulation bounds
*       wL, wU, wtheta, wphi,
*       wBACK                = linearised products of duals and binaries
*       y_null/y_up/y_low(j) = knockout / upregulation / downregulation
*       y_rev(j)             = 1 if the reaction half carries no flux
*****************************************************************************************
Variables
z, b
v(j)
lambda(i)
mu(j)
wL(j), wU(j)
wphi(j), wtheta(j)
wFOR(j), wBACK(j)
;

Negative Variables
biomass
;

Positive variables
phi(j), theta(j)
muL(j), muU(j)
mFor(j),mBack(j)
;

binary variables
y_null(j)
y_low(j)
y_up(j)
y_rev(j)
;

scalar counter /1/;
scalar M /100/;
scalar k ;
scalar n /0/;
scalar cutoff /0/;
scalar cutoff1 /0/;

***************************** EQUATIONS *************************************************
*  Outer level   Obj, Tot, y_lim, Int_cut, Primal_dual (strong duality), big-M
*                linearisations (Lin_muL*, Lin_muU*, UP_outer*, DOWN_outer*, back_lin*)
*  Primal        Stoic, Flux_limL/U, Flux_Upregulate, Flux_Downregulate,
*                RevRXN, ForBackRxn
*  Dual          dual1a, dual1b, dual1c
*****************************************************************************************
Equations

Obj
Stoic
bio, Flux_Upper, Flux_Lower, Flux_Upper_t, Flux_Lower_t
Con_bio
Flux_limL,Flux_limU
Flux_Upregulate
Flux_Downregulate

Primal_dual
dual1a,dual1b,dual1c

Int_cut
Tot
y_lim

back_lin1
back_lin2
back_lin3
back_lin4
RevRXN
ForBackRxn

UP_outer4, UP_outer5, UP_outer6, UP_outer7
DOWN_outer8, DOWN_outer9, DOWN_outer10, DOWN_outer11
Lin_muLa,Lin_muLb,Lin_muLc,Lin_muLd
Lin_muUa,Lin_muUb,Lin_muUc,Lin_muUd
;

bio..                   b =e= v('R3930');
Flux_Upper(j)..         v(j) =l= UB(j);
Flux_Lower(j)..         -v(j) =l= -LB(j);

* bounds with a found regulation set applied, for the FBA check (TestModel)
Flux_Upper_t(j)..       v(j) =l= (1-null_par(j)) * ( UB(j) - (((UB(j) - LB(j))*e_down)*down_par(j)));
Flux_Lower_t(j)..       -v(j) =l= -(1-null_par(j)) * (LB(j) + (((UB(j) - LB(j))*e_up)*up_par(j)));

Con_bio..               -v('R3930') =l= -maxb*cb ;
Obj..                   z =e= v('E3934');
Stoic(i)..              sum(j, (S(i,j)*v(j))) =e= 0;

Flux_Downregulate(j)..  v(j) =l= (UB(j) - ((UB(j) - LB(j))*e_down))*y_low(j) + UB(j)*(1-y_low(j)-y_null(j));
Flux_limL(j)..          -v(j) =l= -LB(j)*(1-y_null(j) );
Flux_limU(j)..          v(j) =l= UB(j)*(1-y_null(j) );
Flux_Upregulate(j)..    -v(j) =l= -((LB(j) + ((UB(j) - LB(j))*e_up))*y_up(j) + LB(j)*(1-y_up(j)-y_null(j)));

Tot..                   sum(j, y_null(j) + y_low(j)+ y_up(j)) =e= k;
y_lim(j)..              y_null(j) + y_low(j) + y_up(j) =l= 1;
Int_cut(index)..        sum(j, done_up(index,j) * y_up(j) + done_ko(index,j) * y_null(j)  + done_down(index,j) * y_low(j)) =l= k-1;

* exactly one half of each reversible reaction is switched off
ForBackRxn(j,jj)$( (Reverse(j) > 0) and (ord(jj)=Reverse(j)) )..  y_rev(j) + y_rev(jj) =e= 1;
RevRXN(j)..             v(j) =l= (1-(y_rev(j))) * UB(j);
back_lin1(j)..          wBACK(j) =l= M*(y_rev(j));
back_lin2(j)..          wBACK(j) =g= -M*(y_rev(j));
back_lin3(j)..          wBACK(j) =l= mBack(j) + M*(1-y_rev(j));
back_lin4(j)..          wBACK(j) =g= mBack(j) - M*(1-y_rev(j));

* strong duality: inner primal objective = inner dual objective
$ifThen %algorithm% == Bio_Igg
Primal_dual..           v('R3930') =e=
$else
Primal_dual..           - v('E3934') =e= maxb*cb*mu('R3930')
$endIf
                        - sum(j, wtheta(j)*(LB(j) + ((UB(j) - LB(j))*e_up)) + theta(j)*LB(j) - wtheta(j)*LB(j) )
                        + sum(j, wphi(j)*(UB(j) - ((UB(j) - LB(j))*e_down)) + phi(j)*UB(j) - wphi(j)*UB(j) )
                        + sum(j, muU(j)*UB(j) - wU(j)*UB(j))
                        - sum(j, muL(j)*LB(j) - wL(j)*LB(j))
                        + sum(j, mBack(j)*UB(j) - wBACK(j)*UB(j));

Lin_muLa(j)..           wL(j) =l= M*y_null(j);
Lin_muLb(j)..           wL(j) =g= -M*y_null(j);
Lin_muLc(j)..           wL(j) =l= muL(j) + M*(1-y_null(j) );
Lin_muLd(j)..           wL(j) =g= muL(j) - M*(1-y_null(j) );

Lin_muUa(j)..           wU(j) =l= M*y_null(j);
Lin_muUb(j)..           wU(j) =g= -M*y_null(j);
Lin_muUc(j)..           wU(j) =l= muU(j) + M*(1-y_null(j) );
Lin_muUd(j)..           wU(j) =g= muU(j) - M*(1-y_null(j) );

DOWN_outer8(j)..        wphi(j) =l= M*y_low(j);
DOWN_outer9(j)..        wphi(j) =g= -M*y_low(j);
DOWN_outer10(j)..       wphi(j) =l= phi(j) + M*(1-y_low(j) );
DOWN_outer11(j)..       wphi(j) =g= phi(j) - M*(1-y_low(j) );

UP_outer4(j)..          wtheta(j) =l= M*y_up(j);
UP_outer5(j)..          wtheta(j) =g= -M*y_up(j);
UP_outer6(j)..          wtheta(j) =l= theta(j) + M*(1-y_up(j) );
UP_outer7(j)..          wtheta(j) =g= theta(j) - M*(1-y_up(j) );

dual1b(j)$(not biom(j))..  sum(i,lambda(i)*S(i,j)) - theta(j) + phi(j) + muU(j) - muL(j) + mBack(j) =e= 0;
dual1c(j)$(biom(j))..      sum(i, lambda(i)*S(i,j)) - theta(j) + phi(j) + muU(j) - muL(j) + mBack(j) =e= 1;
* used by Igg_Igg only
dual1a..                   sum(i, lambda(i)*S(i, 'R3930')) + mu('R3930') + theta('R3930') - phi('R3930') + muU('R3930') - muL('R3930') + mBack('R3930') =e= 0;

************** BOUNDS AND FIXED BINARIES ************************************************

scalar max /100/;

LB(j)$(constraint_flux(j)) = 0;
UB(j)$(constraint_flux(j)) = 10;

v.lo(j)$(not constraint_flux(j) ) = 0;
v.up(j)$(not constraint_flux(j) ) = max;
v.lo(j)$(constraint_flux(j) ) = 0;
v.up(j)$(constraint_flux(j) ) = max;

* irreversible and exchange reactions have no reversibility binary
y_rev.fx(j)$(rxntype(j)=0 or rxntype(j)=3 or rxntype(j)=4) = 0;

* not regulated: medium uptakes, protected reactions, exchanges,
* blocked reactions and reactions without a gene association
y_up.fx(j)$(exchange(j)) = 0;
y_null.fx(j)$(exchange(j)) = 0;
y_low.fx(j)$(exchange(j)) = 0;

y_up.fx(j)$(constraint_flux(j)) = 0;
y_null.fx(j)$(constraint_flux(j)) = 0;
y_low.fx(j)$(constraint_flux(j)) = 0;

y_up.fx(j)$(rxntype(j)=4) = 0;
y_null.fx(j)$(rxntype(j)=4) = 0;
y_low.fx(j)$(rxntype(j)=4) = 0;

y_up.fx(j)$(LB(j)=0 and UB(j)=0) = 0;
y_low.fx(j)$(LB(j)=0 and UB(j)=0) = 0;
y_null.fx(j)$(LB(j)=0 and UB(j)=0) = 0;

y_up.fx(j)$(NoGenes(j)) = 0;
y_low.fx(j)$(NoGenes(j)) = 0;
y_null.fx(j)$(NoGenes(j)) = 0;

* downregulation disabled (see e_down above)
y_low.fx(j) = 0;

********************************** MODELS ***********************************************

Model MaxBio
/
Stoic bio Flux_Upper Flux_Lower
/;

Model MaxIgG
/
Stoic Obj Flux_Upper Flux_Lower
/;

* FBA check of a found regulation set
Model TestModel
/
Stoic bio Flux_Upper_t Flux_Lower_t
/

model bilevel
/
Obj
Stoic
$ifThen %algorithm% == Igg_Igg
Con_bio
dual1a
$endIf

Flux_limL,Flux_limU
Flux_Upregulate
Flux_Downregulate

Primal_dual
dual1b
dual1c

Int_cut
Tot
y_lim

back_lin1
back_lin2
back_lin3
back_lin4
RevRXN
ForBackRxn

UP_outer4, UP_outer5, UP_outer6, UP_outer7
DOWN_outer8, DOWN_outer9, DOWN_outer10, DOWN_outer11
Lin_muLa,Lin_muLb,Lin_muLc,Lin_muLd
Lin_muUa,Lin_muUb,Lin_muUc,Lin_muUd
/
;

*************** SOLVE LOOP **************************************************************
*  For each (e_up, k): maximise IgG, record the regulation set, add an integer cut
*  and solve again, up to five alternative solutions. A solution is written out only
*  if the FBA check (TestModel) with that regulation set still gives IgG > 0.
*****************************************************************************************

file forced /RXNOptimization_%doc%_%algorithm%.txt/;
put forced;
put '****************************************'//;

solve MaxBio using lp maximizing b;
solve MaxIgG using lp maximizing z;
put 'Max Biomass ='b.l:0:8/;
put 'Max IgG ='z.l:0:8/;
maxb = b.l;
maxi = z.l;

* biomass must stay at or above 10% of its maximum
$ifThen %algorithm% == Bio_Igg
LB(j)$(biom(j)) = %maxgrowth% * 0.1;
$endIf
cb = 0.1;
ci = 0;

e_down = %edown%;

for (e_up = %eup_from% to %eup_to% by %eup_step%,

cutoff = 0;
cutoff1 = 0;
jlnj = 0;

for (k = %k_from% to %k_to% by %k_step%,

put "Upregulation coefficient is: "e_up, "Downregulation coefficient is: ",e_down /;

    counter = 1;
    n = 0;

    done_up(index,j) = 0;
    done_ko(index,j) = 0;
    done_down(index,j) = 0;

    put 'Number of Interventions (k) =' k/;
    put '**********************************************************'//;
    while( counter = 1,

        solve bilevel using mip maximizing z;
        jlnj = jlnj + 1;
        n = n + 1;
        if(n eq 1, cutoff = z.l );

        if(bilevel.modelstat eq 1 or bilevel.modelstat eq 8 ,
            if( z.l <= cutoff1  * 1.01,
                counter = 0;
                put "Objective value/better solution not found Cutoff1 = ", cutoff1:0:8, "z = ", z.l:0:8 /;
                cutoff1 = cutoff;
            );

            put 'Total Iterations are' jlnj //;
            put 'Biomass = 'v.l('R3930'):0:8/;
            put 'Antibody = 'v.l('E3934'):0:8/;

            up_par(j) = y_up.l(j);
            null_par(j) = y_null.l(j);
            down_par(j) = y_low.l(j);
            solve TestModel using lp maximizing b;

            if ( v.l('E3934') > 0,
                put 'Successful Iteration No. =' n//;
                put 'Biomass = 'v.l('R3930'):0:8/;
                put 'Antibody = 'v.l('E3934'):0:8/;
                put '****************************************'//;
                put /'Gene Regulations:'//;
                put "%%%%% Increase the:"/;
                w1 = 0;
                loop (j$(y_up.l(j) = 1),
                    w1 = w1 + 1;
                    put "r",w1:0:0,"=",ord(j):0:0,";"/;
                );
                put "b = [];"/"for i=1:"w1:0:0/"eval(sprintf('b = [b; r%d];',i));"/"end"/"r_up = b;"/;
                w2 = w1;
                w1 = w1 + 1;
                put "%%%%% Knock-out the:"/;
                loop (j$(y_null.l(j) = 1),
                    w2 = w2 + 1;
                    put "r",w2:0:0,"=",ord(j):0:0,";"/;
                );
                put "b = [];"/"for i="w1:0:0":"w2:0:0/"eval(sprintf('b = [b; r%d];',i));"/"end"/"r_null = b;"/;
            );

            done_up(index,j)$(ord(index) = n and (y_up.l(j) = 1)) = 1;
            done_ko(index,j)$(ord(index) = n and (y_null.l(j) = 1)) = 1;
            done_down(index,j)$(ord(index) = n and (y_low.l(j) = 1)) = 1;
        );

        if(n >= 5 or (bilevel.modelstat ne 8 and bilevel.modelstat ne 1),
            counter = 0;
            put "No feasible solution achieved, modelstat :"bilevel.modelstat, " n = ",n//;
            cutoff1 = cutoff;
        );

    put '**************************************************'//;
    );

    put /'**********************************************************'/;
);
);

putclose forced;
*****************************************************************************************
