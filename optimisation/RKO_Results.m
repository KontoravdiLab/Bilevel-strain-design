 function [up,r_up,r_null,Algorithm,Feed,Regulations] =  RKO_Results(x,model,model_type)
 
%  if model_type == 'Secretome'
%      model_type = 'K1';
%  end

if (strcmp(model_type,'K1'))
    if (x == 1)
        
      Regulations = 50;
      Feed = 'D1';
      up = 0.5;
      Algorithm = 'Bio_Igg';
      r_up(1) = find(ismember(model.rxns,'Htx_f'));
      r_up(2) = find(ismember(model.rxns,'PPItx_f'));
      r_up(3) = find(ismember(model.rxns,'PIt2m_2'));
      r_up(4) = find(ismember(model.rxns,'NADtx_b'));
      r_null(1) = find(ismember(model.rxns,'GAPD_f'));
      r_null(2) = find(ismember(model.rxns,'PPAm'));
      r_null(3) = find(ismember(model.rxns,'34HPPte_r'));
      r_null(4) = find(ismember(model.rxns,'Htx_b'));
  elseif (x == 2)
      Regulations = 60;
      Feed = 'D1';
      up = 0.4;
      Algorithm = 'Bio_Igg';
      r_up(1) = find(ismember(model.rxns,'PItn_f'));
      r_up(2) = find(ismember(model.rxns,'PPItx_f'));
      r_up(3) = find(ismember(model.rxns,'GLYC3Ptmc_b'));
      r_up(4) = find(ismember(model.rxns,'H2Otly_b'));
      r_up(5) = find(ismember(model.rxns,'r1290_b'));
      r_up(6) = find(ismember(model.rxns,'RTOTAL3t_b'));
      r_null(1) = find(ismember(model.rxns,'34HPPte_r'));
      r_null(2) = find(ismember(model.rxns,'GLYC3Ptmc_f'));
      r_null(3) = find(ismember(model.rxns,'H2Otly_f'));
      r_null(4) = find(ismember(model.rxns,'r1290_f'));
  elseif (x == 3)
      Regulations = 50;
      Feed = 'D1';
      up = 0.6;
      Algorithm = 'Bio_Igg';
      r_up(1) = find(ismember(model.rxns,'PPAm'));
      r_null(1) = find(ismember(model.rxns,'GLXO1'));
      r_null(2) = find(ismember(model.rxns,'GUAt_r'));
      r_null(3) = find(ismember(model.rxns,'34HPPte_r'));
  elseif (x == 4)
      Regulations = 40;
      Feed = 'D4';
      up = 0.5;
      Algorithm = 'Igg_Igg';
      r_up(1) = find(ismember(model.rxns,'DESAT16_2'));
      r_up(2) = find(ismember(model.rxns,'DGAT_cho'));
      r_up(3) = find(ismember(model.rxns,'CDS_cho'));
      r_up(4) = find(ismember(model.rxns,'ACGAM6PSi'));
      r_up(5) = find(ismember(model.rxns,'PETHCT'));
      r_null(1) = find(ismember(model.rxns,'RE3273C'));
      r_null(2) = find(ismember(model.rxns,'PROD2'));
      r_null(3) = find(ismember(model.rxns,'PIt2m_f'));
      r_null(4) = find(ismember(model.rxns,'ASPTAm_f'));
      r_null(5) = find(ismember(model.rxns,'NADPtxu'));
      r_null(6) = find(ismember(model.rxns,'peplys_synthesis'));
      r_null(7) = find(ismember(model.rxns,'SUCOAS1m_b'));
      r_null(8) = find(ismember(model.rxns,'SUCOASm_b'));
  elseif (x == 5)
      Regulations = 50;
      Feed = 'D4';
      up = 0.5;
      Algorithm = 'Igg_Igg';
      r_up(1) = find(ismember(model.rxns,'DGAT_cho'));
      r_up(2) = find(ismember(model.rxns,'3DSPHR'));
      r_up(3) = find(ismember(model.rxns,'O2ter'));
      r_null(1) = find(ismember(model.rxns,'THRA'));
      r_null(2) = find(ismember(model.rxns,'PYRtm_f'));
      r_null(3) = find(ismember(model.rxns,'BTNPL'));
      r_null(4) = find(ismember(model.rxns,'PEPCKm'));
      r_null(5) = find(ismember(model.rxns,'LYStn_f'));
      r_null(6) = find(ismember(model.rxns,'PPItm_f'));
  elseif (x == 6)
      Regulations = 40;
      Feed = 'D4';
      up = 0.7;
      Algorithm = 'Igg_Igg';
      r_up(1) = find(ismember(model.rxns,'GLNS'));
      r_up(2) = find(ismember(model.rxns,'GPAM_cho'));
      r_up(3) = find(ismember(model.rxns,'URIDK2_f'));
      r_up(4) = find(ismember(model.rxns,'PPM_f'));
      r_up(5) = find(ismember(model.rxns,'PETHCT'));
      r_null(1) = find(ismember(model.rxns,'PIt2m_f'));
      r_null(2) = find(ismember(model.rxns,'VALTA'));
      r_null(3) = find(ismember(model.rxns,'CITRt2m_b'));
  elseif (x == 7)
      Regulations = 45;
      Feed = 'D4';
      up = 0.8;
      Algorithm = 'Igg_Igg';
      r_up(1) = find(ismember(model.rxns,'PSSA1_cho_f'));
      r_up(2) = find(ismember(model.rxns,'FUM_f'));
      r_up(3) = find(ismember(model.rxns,'NDP7g'));
      r_up(4) = find(ismember(model.rxns,'O2ter'));
      r_null(1) = find(ismember(model.rxns,'ACACT1rm_f'));
      r_null(2) = find(ismember(model.rxns,'34HPPte_r'));
      r_null(3) = find(ismember(model.rxns,'GLCtly'));
      r_null(4) = find(ismember(model.rxns,'NADPHtxu'));
      r_null(5) = find(ismember(model.rxns,'PPItm_f'));
      r_null(6) = find(ismember(model.rxns,'ASPTA_b'));
      r_null(7) = find(ismember(model.rxns,'RPE_b'));
  elseif (x == 8)
      Regulations = 40;
      Feed = 'D4';
      up = 0.9;
      Algorithm = 'Igg_Igg';
      r_up(1) = find(ismember(model.rxns,'DGAT_cho'));
      r_up(2) = find(ismember(model.rxns,'URIDK2_f'));
      r_up(3) = find(ismember(model.rxns,'GF6PTA'));
      r_null(1) = find(ismember(model.rxns,'PROD2'));
      r_null(2) = find(ismember(model.rxns,'THRA'));
      r_null(3) = find(ismember(model.rxns,'CITRtm_f'));
      r_null(4) = find(ismember(model.rxns,'THRD_L'));
      r_null(5) = find(ismember(model.rxns,'PPA'));
      r_null(6) = find(ismember(model.rxns,'DUTPDPm'));
      r_null(7) = find(ismember(model.rxns,'PPItm_f'));
      r_null(8) = find(ismember(model.rxns,'PPItx_f'));
      r_null(9) = find(ismember(model.rxns,'THR3D'));
      r_null(10) = find(ismember(model.rxns,'ASPGLUm_b'));
  elseif (x == 9)
      Regulations = 35;
      Feed = 'D1';
      up = 0.9;
      Algorithm = 'Igg_Igg';
      r_up(1) = find(ismember(model.rxns,'URIDK2_f'));
      r_up(2) = find(ismember(model.rxns,'PPM_f'));
      r_up(3) = find(ismember(model.rxns,'IMPD'));
      r_up(4) = find(ismember(model.rxns,'FUMm_b'));
      r_null(1) = find(ismember(model.rxns,'ACACT1rm_f'));
      r_null(2) = find(ismember(model.rxns,'34HPPte_r'));
      r_null(3) = find(ismember(model.rxns,'NADPtxu'));
      r_null(4) = find(ismember(model.rxns,'PPItm_f'));
      r_null(5) = find(ismember(model.rxns,'PROtm_f'));
      r_null(6) = find(ismember(model.rxns,'r0081_f'));
      r_null(7) = find(ismember(model.rxns,'r0202_b'));
  elseif (x == 10)
      Regulations = 60;
      Feed = 'D1';
      up = 0.9;
      Algorithm = 'Igg_Igg';
      r_up(1) = find(ismember(model.rxns,'DESAT16_2'));
      r_up(2) = find(ismember(model.rxns,'URIDK2_f'));
      r_up(3) = find(ismember(model.rxns,'PPM_f'));
      r_up(4) = find(ismember(model.rxns,'CDIPTr_cho'));
      r_up(5) = find(ismember(model.rxns,'IMPC'));
      r_up(6) = find(ismember(model.rxns,'DGK1_b'));
      r_null(1) = find(ismember(model.rxns,'IZPN'));
      r_null(2) = find(ismember(model.rxns,'NDP8'));
      r_null(3) = find(ismember(model.rxns,'PPAm'));
      r_null(4) = find(ismember(model.rxns,'34HPPte_r'));
      r_null(5) = find(ismember(model.rxns,'DOLP_ter'));
      r_null(6) = find(ismember(model.rxns,'PROtm_f'));
      r_null(7) = find(ismember(model.rxns,'L_LACtm_b'));
      r_null(8) = find(ismember(model.rxns,'SUCOAS1m_b'));
      r_null(9) = find(ismember(model.rxns,'SUCOASm_b'));
      r_null(10) = find(ismember(model.rxns,'r0424_b'));
      r_null(11) = find(ismember(model.rxns,'r0480_b'));
  elseif (x == 11)
      Regulations = 35;
      Feed = 'D1';
      up = 1;
      Algorithm = 'Igg_Igg';
      r_up(1) = find(ismember(model.rxns,'GPAM_cho'));
      r_null(1) = find(ismember(model.rxns,'34HPPte_r'));
  elseif (x == 12)
      Regulations = 50;
      Feed = 'D1';
      up = 1;
      Algorithm = 'Igg_Igg';
      r_up(1) = find(ismember(model.rxns,'PPM_f'));
      r_null(1) = find(ismember(model.rxns,'SMOX'));
      r_null(2) = find(ismember(model.rxns,'PIt2m_f'));
      r_null(3) = find(ismember(model.rxns,'GAPD_f'));
      r_null(4) = find(ismember(model.rxns,'34HPPte_r'));
  elseif (x == 13)
      Regulations = 35;
      Feed = 'D4';
      up = 0.8;
      Algorithm = 'Bio_Igg';
      r_up(1) = find(ismember(model.rxns,'AHCYSte_r'));
      r_null(1) = find(ismember(model.rxns,'PEPLYStn_r'));
  elseif (x == 14)
      Regulations = 35;
      Feed = 'D4';
      up = 0.9;
      Algorithm = 'Bio_Igg';
      r_up(1) = find(ismember(model.rxns,'ACGAMPM'));
      r_up(2) = find(ismember(model.rxns,'URIt_f'));
      r_up(3) = find(ismember(model.rxns,'GLU5SAtmc_f'));
      r_null(1) = find(ismember(model.rxns,'GLU5SAtmc_b'));
      r_null(2) = find(ismember(model.rxns,'PPItr_b'));
      r_null(3) = find(ismember(model.rxns,'r2521_b'));
  elseif (x == 15)
      Regulations = 50;
      Feed = 'D3';
      up = 1;
      Algorithm = 'Igg_Igg';
      r_null(1) = find(ismember(model.rxns,'PROD2'));
      r_null(2) = find(ismember(model.rxns,'THRA'));
      r_null(3) = find(ismember(model.rxns,'PIt2m_f'));
      r_null(4) = find(ismember(model.rxns,'LEUTA'));
      r_null(5) = find(ismember(model.rxns,'VALTA'));
      r_null(6) = find(ismember(model.rxns,'34HPPte_r'));
      r_null(7) = find(ismember(model.rxns,'ASNtm'));
      r_null(8) = find(ismember(model.rxns,'FADtru'));
      r_null(9) = find(ismember(model.rxns,'ACACt2m_b'));
      r_null(10) = find(ismember(model.rxns,'PYRt2m_b'));
      r_null(11) = find(ismember(model.rxns,'r0911_b'));
      r_null(12) = find(ismember(model.rxns,'FADRx_b'));
  elseif (x == 16)
      Regulations = 60;
      Feed = 'D3';
      up = 1;
      Algorithm = 'Igg_Igg';
      r_null(1) = find(ismember(model.rxns,'PYRtm_f'));
      r_null(2) = find(ismember(model.rxns,'ACACT1rm_f'));
      r_null(3) = find(ismember(model.rxns,'ASPTAm_f'));
      r_null(4) = find(ismember(model.rxns,'LEUTA'));
      r_null(5) = find(ismember(model.rxns,'34HPPte_r'));
      r_null(6) = find(ismember(model.rxns,'ARTPLM1_f'));
      r_null(7) = find(ismember(model.rxns,'CO2tm_f'));
      r_null(8) = find(ismember(model.rxns,'O2tm_f'));
      r_null(9) = find(ismember(model.rxns,'r0947_f'));
  elseif (x == 17)
      Regulations = 60;
      Feed = 'D1';
      up = 0.5;
      Algorithm = 'Bio_Igg';
      r_up(1) = find(ismember(model.rxns,'CAT'));
      r_up(2) = find(ismember(model.rxns,'MDHm_f'));
      r_up(3) = find(ismember(model.rxns,'ME2'));
      r_up(4) = find(ismember(model.rxns,'P5CDm'));
      r_up(5) = find(ismember(model.rxns,'PItx_r'));
      r_up(6) = find(ismember(model.rxns,'r2384'));
      r_null(1) = find(ismember(model.rxns,'34HPPte_r'));
  elseif (x == 18)
      Regulations = 47;
      Feed = 'D4';
      up = 0.7;
      Algorithm = 'Bio_Igg';
      r_up(1) = find(ismember(model.rxns,'PPAm'));
      r_up(2) = find(ismember(model.rxns,'r1290_f'));
      r_null(1) = find(ismember(model.rxns,'FADRx_b'));
  elseif (x == 19)
      Regulations = 30;
      Feed = 'D4';
      up = 0.7;
      Algorithm = 'Bio_Igg';
      r_up(1) = find(ismember(model.rxns,'LEUTA'));
      r_null(1) = find(ismember(model.rxns,'G5SADr_b'));
  elseif (x == 20)
      Regulations = 200;
      Feed = 'D2';
      up = 0.7;
      Algorithm = 'Bio_Igg';
      r_up(1) = find(ismember(model.rxns,'CAT'));
      r_up(2) = find(ismember(model.rxns,'r0383'));
      r_up(3) = find(ismember(model.rxns,'r1303'));
      r_up(4) = find(ismember(model.rxns,'r2465'));
      r_null(1) = find(ismember(model.rxns,'LEUTA'));
      r_null(2) = find(ismember(model.rxns,'KHte_f'));
      r_null(3) = find(ismember(model.rxns,'NAt3_1_b'));
      r_null(4) = find(ismember(model.rxns,'OCOAT1m_b'));
      r_null(5) = find(ismember(model.rxns,'PPItr_b'));
      r_null(6) = find(ismember(model.rxns,'r2521_b'));
  elseif (x == 21)
      Regulations = 90;
      Feed = 'D4';
      up = 0.7;
      Algorithm = 'Bio_Igg';
      r_up(1) = find(ismember(model.rxns,'NADtru'));
      r_up(2) = find(ismember(model.rxns,'r0062'));
      r_up(3) = find(ismember(model.rxns,'PIt2m_b'));
      r_null(1) = find(ismember(model.rxns,'PROD2'));
      r_null(2) = find(ismember(model.rxns,'PIt2m_f'));
      r_null(3) = find(ismember(model.rxns,'r0202_b'));
      r_null(4) = find(ismember(model.rxns,'FADRx_b'));
  elseif (x == 22)
      Regulations = 90;
      Feed = 'D4';
      up = 0.7;
      Algorithm = 'Bio_Igg';
      r_up(1) = find(ismember(model.rxns,'FUM_f'));
      r_up(2) = find(ismember(model.rxns,'H2Otp_f'));
      r_up(3) = find(ismember(model.rxns,'Htx_f'));
      r_up(4) = find(ismember(model.rxns,'r0062'));
      r_up(5) = find(ismember(model.rxns,'r1290_f'));
      r_null(1) = find(ismember(model.rxns,'PROD2'));
      r_null(2) = find(ismember(model.rxns,'PYRtm_f'));
      r_null(3) = find(ismember(model.rxns,'Htx_b'));
      r_null(4) = find(ismember(model.rxns,'r0202_b'));
  elseif (x == 23)
      Regulations = 90;
      Feed = 'D4';
      up = 0.8;
      Algorithm = 'Bio_Igg';
      r_up(1) = find(ismember(model.rxns,'r0165'));
      r_up(2) = find(ismember(model.rxns,'r1502'));
      r_up(3) = find(ismember(model.rxns,'TPI_b'));
      r_up(4) = find(ismember(model.rxns,'RTOTAL3t_b'));
      r_null(1) = find(ismember(model.rxns,'PROD2'));
  elseif (x == 24)
      Regulations = 80;
      Feed = 'D4';
      up = 0.9;
      Algorithm = 'Bio_Igg';
      r_up(1) = find(ismember(model.rxns,'PPAm'));
      r_null(1) = find(ismember(model.rxns,'FAS140ACP'));
      r_null(2) = find(ismember(model.rxns,'KAS8'));
      r_null(3) = find(ismember(model.rxns,'CHOLtu_f'));
      r_null(4) = find(ismember(model.rxns,'PROD2'));
      r_null(5) = find(ismember(model.rxns,'PPItr_b'));
      r_null(6) = find(ismember(model.rxns,'r2521_b'));
  elseif (x == 25)
      Regulations = 90;
      Feed = 'D4';
      up = 0.9;
      Algorithm = 'Bio_Igg';
      r_up(1) = find(ismember(model.rxns,'PPItm_f'));
      r_up(2) = find(ismember(model.rxns,'r1502'));
      r_null(1) = find(ismember(model.rxns,'PROD2'));
      r_null(2) = find(ismember(model.rxns,'FUM_b'));
  elseif (x == 26)
      Regulations = 70;
      Feed = 'D4';
      up = 1.1;
      Algorithm = 'Bio_Igg';
      r_null(1) = find(ismember(model.rxns,'HMGLm'));
      r_null(2) = find(ismember(model.rxns,'PIt2m_f'));
      r_null(3) = find(ismember(model.rxns,'CSm'));
      r_null(4) = find(ismember(model.rxns,'GLPASE2a'));
      r_null(5) = find(ismember(model.rxns,'ACOATA'));
      r_null(6) = find(ismember(model.rxns,'MI14PP'));
      r_null(7) = find(ismember(model.rxns,'VALTA'));
      r_null(8) = find(ismember(model.rxns,'r0047'));
      r_null(9) = find(ismember(model.rxns,'peplys_synthesis'));
      r_null(10) = find(ismember(model.rxns,'ASPTA_b'));
      r_null(11) = find(ismember(model.rxns,'FUM_b'));
      r_null(12) = find(ismember(model.rxns,'TALA_b'));
      r_null(13) = find(ismember(model.rxns,'r1109_b'));
  elseif (x == 27)
      Regulations = 70;
      Feed = 'D1';
      up = 0.7;
      Algorithm = 'Bio_Igg';
      r_up(1) = find(ismember(model.rxns,'CLHCOtex'));
      r_up(2) = find(ismember(model.rxns,'PIt2m_f'));
      r_up(3) = find(ismember(model.rxns,'FADH2ETC'));
      r_up(4) = find(ismember(model.rxns,'r1492'));
      r_null(1) = find(ismember(model.rxns,'AACOAT'));
      r_null(2) = find(ismember(model.rxns,'34HPPte_r'));
      r_null(3) = find(ismember(model.rxns,'KHte_f'));
      r_null(4) = find(ismember(model.rxns,'r0082_b'));
      r_null(5) = find(ismember(model.rxns,'r2137_b'));
  elseif (x == 28)
      Regulations = 70;
      Feed = 'D3';
      up = 1.1;
      Algorithm = 'Igg_Igg';
      r_null(1) = find(ismember(model.rxns,'SMOX'));
      r_null(2) = find(ismember(model.rxns,'FBA'));
      r_null(3) = find(ismember(model.rxns,'HMGCOASim'));
      r_null(4) = find(ismember(model.rxns,'LEUTA'));
      r_null(5) = find(ismember(model.rxns,'PGCD'));
      r_null(6) = find(ismember(model.rxns,'SBTD_D2'));
      r_null(7) = find(ismember(model.rxns,'34HPPte_r'));
      r_null(8) = find(ismember(model.rxns,'4HPRO_LTte'));
      r_null(9) = find(ismember(model.rxns,'ARTPLM1_f'));
      r_null(10) = find(ismember(model.rxns,'r0407'));
      r_null(11) = find(ismember(model.rxns,'GLU5SAtmc_b'));
      r_null(12) = find(ismember(model.rxns,'gthox_export_b'));
  elseif (x == 29)
      Regulations = 100;
      Feed = 'D2';
      up = 0.7;
      Algorithm = 'Igg_Igg';
      r_up(1) = find(ismember(model.rxns,'ACACt2m_f'));
      r_up(2) = find(ismember(model.rxns,'r0082_f'));
      r_up(3) = find(ismember(model.rxns,'r1503'));
      r_up(4) = find(ismember(model.rxns,'TPI_b'));
      r_up(5) = find(ismember(model.rxns,'H2Otn_b'));
      r_up(6) = find(ismember(model.rxns,'r2179_b'));
      r_up(7) = find(ismember(model.rxns,'r2180_b'));
      r_up(8) = find(ismember(model.rxns,'RTOTAL3t_b'));
      r_null(1) = find(ismember(model.rxns,'VALTA'));
      r_null(2) = find(ismember(model.rxns,'NADtx_f'));
      r_null(3) = find(ismember(model.rxns,'r1290_f'));
      r_null(4) = find(ismember(model.rxns,'r2178_f'));
      r_null(5) = find(ismember(model.rxns,'r2179_f'));
      r_null(6) = find(ismember(model.rxns,'MDHm_b'));
      r_null(7) = find(ismember(model.rxns,'LDH_Lm_b'));
      r_null(8) = find(ismember(model.rxns,'r2516_b'));
  elseif (x == 30)
      Regulations = 110;
      Feed = 'D2';
      up = 0.8;
      Algorithm = 'Igg_Igg';
      r_up(1) = find(ismember(model.rxns,'SUCOASm_f'));
      r_up(2) = find(ismember(model.rxns,'GLU5SAtmc_f'));
      r_up(3) = find(ismember(model.rxns,'Htx_f'));
      r_up(4) = find(ismember(model.rxns,'RTOT1_cho'));
      r_up(5) = find(ismember(model.rxns,'r2516_b'));
      r_null(1) = find(ismember(model.rxns,'GLXO1'));
      r_null(2) = find(ismember(model.rxns,'PIK4_cho'));
      r_null(3) = find(ismember(model.rxns,'G3PD2m'));
      r_null(4) = find(ismember(model.rxns,'PGCD'));
      r_null(5) = find(ismember(model.rxns,'DHAPtc'));
      r_null(6) = find(ismember(model.rxns,'r2516_f'));
      r_null(7) = find(ismember(model.rxns,'HMGCOARc_b'));
      r_null(8) = find(ismember(model.rxns,'GLU5SAtmc_b'));
  elseif (x == 31)
      Regulations = 110;
      Feed = 'D2';
      up = 0.9;
      Algorithm = 'Igg_Igg';
      r_up(1) = find(ismember(model.rxns,'G5SADr_f'));
      r_null(1) = find(ismember(model.rxns,'G5SADr_b'));
  elseif (x == 32)
      Regulations = 200;
      Feed = 'D3';
      up = 0.9;
      Algorithm = 'Bio_Igg';
      r_up(1) = find(ismember(model.rxns,'G5SADr_f'));
      r_null(1) = find(ismember(model.rxns,'G5SADr_b'));
  elseif (x == 33)
      Regulations = 90;
      Feed = 'D3';
      up = 0.9;
      Algorithm = 'Bio_Igg';
      r_up(1) = find(ismember(model.rxns,'GLU5SAtmc_f'));
      r_null(1) = find(ismember(model.rxns,'GLU5SAtmc_b'));
  elseif (x == 34)
      Regulations = 80;
      Feed = 'D2';
      up = 0.8;
      Algorithm = 'Bio_Igg';
      r_up(1) = find(ismember(model.rxns,'UGCGc'));
      r_up(2) = find(ismember(model.rxns,'SUCOASm_f'));
      r_up(3) = find(ismember(model.rxns,'r1531'));
      r_null(1) = find(ismember(model.rxns,'PGCD'));
      r_null(2) = find(ismember(model.rxns,'VALTA'));
      r_null(3) = find(ismember(model.rxns,'NADPHtxu'));
      r_null(4) = find(ismember(model.rxns,'r2539'));
      r_null(5) = find(ismember(model.rxns,'ARTPLM2_b'));
  elseif (x == 35)
      Regulations = 50;
      Feed = 'D4';
      up = 0.60;
      Algorithm = 'Igg_Igg';
      r_up(1) = find(ismember(model.rxns,'SPODMm'));
      r_up(2) = find(ismember(model.rxns,'PItn_f'));
      r_up(3) = find(ismember(model.rxns,'r0407'));
      r_up(4) = find(ismember(model.rxns,'r1025'));
      r_up(5) = find(ismember(model.rxns,'r1162'));
      r_up(6) = find(ismember(model.rxns,'r1290_f'));
      r_null(1) = find(ismember(model.rxns,'G5SADr_b'));
      r_null(2) = find(ismember(model.rxns,'r1290_b'));
  elseif (x == 36)
      Regulations = 175;
      Feed = 'D4';
      up = 0.6;
      Algorithm = 'Igg_Igg';
      r_up(1) = find(ismember(model.rxns,'BPNT'));
      r_up(2) = find(ismember(model.rxns,'CSNATm'));
      r_up(3) = find(ismember(model.rxns,'C16txc'));
      r_up(4) = find(ismember(model.rxns,'GLCMter'));
      r_up(5) = find(ismember(model.rxns,'GMPS'));
      r_up(6) = find(ismember(model.rxns,'r0321'));
      r_up(7) = find(ismember(model.rxns,'UDPGNP'));
      r_up(8) = find(ismember(model.rxns,'r1292_b'));
      r_null(1) = find(ismember(model.rxns,'FAS140ACP'));
      r_null(2) = find(ismember(model.rxns,'THRA'));
      r_null(3) = find(ismember(model.rxns,'CSNATp_f'));
      r_null(4) = find(ismember(model.rxns,'PETHCT'));
      r_null(5) = find(ismember(model.rxns,'GLCter_f'));
      r_null(6) = find(ismember(model.rxns,'r1290_f'));
      r_null(7) = find(ismember(model.rxns,'peplys_synthesis'));
  elseif (x == 37)
      Regulations = 225;
      Feed = 'D4';
      up = 0.6;
      Algorithm = 'Igg_Igg';
      r_up(1) = find(ismember(model.rxns,'PYK'));
      r_up(2) = find(ismember(model.rxns,'PRPPS_f'));
      r_up(3) = find(ismember(model.rxns,'PItx_r'));
      r_up(4) = find(ismember(model.rxns,'RTOTAL2t_b'));
      r_null(1) = find(ismember(model.rxns,'PROD2'));
      r_null(2) = find(ismember(model.rxns,'TPI_f'));
      r_null(3) = find(ismember(model.rxns,'COAtm_f'));
      r_null(4) = find(ismember(model.rxns,'PYRt2m_b'));
  elseif (x == 38)
      Regulations = 75;
      Feed = 'D4';
      up = 0.7;
      Algorithm = 'Igg_Igg';
      r_up(1) = find(ismember(model.rxns,'LPS2_cho'));
      r_up(2) = find(ismember(model.rxns,'LACLt'));
      r_up(3) = find(ismember(model.rxns,'PItx_r'));
      r_up(4) = find(ismember(model.rxns,'r0173_f'));
      r_up(5) = find(ismember(model.rxns,'r1290_f'));
      r_null(1) = find(ismember(model.rxns,'ACACT1rm_f'));
      r_null(2) = find(ismember(model.rxns,'Htx_b'));
      r_null(3) = find(ismember(model.rxns,'FADRx_b'));
  elseif (x == 39)
      Regulations = 225;
      Feed = 'D4';
      up = 0.7;
      Algorithm = 'Igg_Igg';
      r_up(1) = find(ismember(model.rxns,'PGK_f'));
      r_up(2) = find(ismember(model.rxns,'HPYRRy'));
      r_up(3) = find(ismember(model.rxns,'SUCOASm_f'));
      r_up(4) = find(ismember(model.rxns,'COAtm_b'));
      r_up(5) = find(ismember(model.rxns,'NADtx_b'));
      r_null(1) = find(ismember(model.rxns,'PROD2'));
      r_null(2) = find(ismember(model.rxns,'r0596'));
      r_null(3) = find(ismember(model.rxns,'LDH_Lm_f'));
      r_null(4) = find(ismember(model.rxns,'COAtm_f'));
      r_null(5) = find(ismember(model.rxns,'ARTPLM1_f'));
      r_null(6) = find(ismember(model.rxns,'NADtx_f'));
      r_null(7) = find(ismember(model.rxns,'r0389_f'));
      r_null(8) = find(ismember(model.rxns,'r1434'));
      r_null(9) = find(ismember(model.rxns,'GLYtm_b'));
  elseif (x == 40)
      Regulations = 125;
      Feed = 'D4';
      up = 0.8;
      Algorithm = 'Igg_Igg';
      r_up(1) = find(ismember(model.rxns,'r0604'));
      r_null(1) = find(ismember(model.rxns,'PROD2'));
  elseif (x == 41)
      Regulations = 200;
      Feed = 'D4';
      up = 0.8;
      Algorithm = 'Igg_Igg';
      r_up(1) = find(ismember(model.rxns,'NKCCt_f'));
      r_up(2) = find(ismember(model.rxns,'G5SADr_f'));
      r_up(3) = find(ismember(model.rxns,'PCRNtm'));
      r_up(4) = find(ismember(model.rxns,'r0160'));
      r_up(5) = find(ismember(model.rxns,'FOLTle'));
      r_up(6) = find(ismember(model.rxns,'H2Otp_f'));
      r_up(7) = find(ismember(model.rxns,'PRISTANALtx'));
      r_up(8) = find(ismember(model.rxns,'r0062'));
      r_up(9) = find(ismember(model.rxns,'r1492'));
      r_up(10) = find(ismember(model.rxns,'KHte_b'));
      r_null(1) = find(ismember(model.rxns,'NAHCO3_HCLt2_f'));
      r_null(2) = find(ismember(model.rxns,'COAtm_f'));
      r_null(3) = find(ismember(model.rxns,'PETHCT'));
      r_null(4) = find(ismember(model.rxns,'OBDHc'));
      r_null(5) = find(ismember(model.rxns,'PItx_r'));
      r_null(6) = find(ismember(model.rxns,'NKCCt_b'));
      r_null(7) = find(ismember(model.rxns,'G5SADr_b'));
      r_null(8) = find(ismember(model.rxns,'CSNATp_b'));
  elseif (x == 42)
      Regulations = 200;
      Feed = 'D4';
      up = 0.8;
      Algorithm = 'Igg_Igg';
      r_up(1) = find(ismember(model.rxns,'CLOXAtex2'));
      r_up(2) = find(ismember(model.rxns,'PSP_L'));
      r_up(3) = find(ismember(model.rxns,'r1487'));
      r_up(4) = find(ismember(model.rxns,'PPItm_f'));
      r_up(5) = find(ismember(model.rxns,'r1501'));
      r_up(6) = find(ismember(model.rxns,'r2465'));
      r_up(7) = find(ismember(model.rxns,'H2Otly_b'));
      r_up(8) = find(ismember(model.rxns,'KHte_b'));
      r_null(1) = find(ismember(model.rxns,'FATP3t_f'));
      r_null(2) = find(ismember(model.rxns,'TRIOK'));
      r_null(3) = find(ismember(model.rxns,'PROD2'));
      r_null(4) = find(ismember(model.rxns,'CHOLK'));
      r_null(5) = find(ismember(model.rxns,'ASNNm'));
      r_null(6) = find(ismember(model.rxns,'ARTPLM1_f'));
      r_null(7) = find(ismember(model.rxns,'KHte_f'));
      r_null(8) = find(ismember(model.rxns,'NADPtxu'));
      r_null(9) = find(ismember(model.rxns,'MMMm_b'));
      r_null(10) = find(ismember(model.rxns,'PE_CHOtm_b'));
  elseif (x == 43)
      Regulations = 200;
      Feed = 'D4';
      up = 0.9;
      Algorithm = 'Igg_Igg';
      r_up(1) = find(ismember(model.rxns,'DURIK1'));
      r_up(2) = find(ismember(model.rxns,'ACACT1rm_f'));
      r_up(3) = find(ismember(model.rxns,'C6DCCACT_r'));
      r_up(4) = find(ismember(model.rxns,'NADPtxu'));
      r_up(5) = find(ismember(model.rxns,'r0156'));
      r_null(1) = find(ismember(model.rxns,'CLOXAtex2'));
      r_null(2) = find(ismember(model.rxns,'G5SADr_b'));
  elseif (x == 44)
      Regulations = 225;
      Feed = 'D4';
      up = 0.9;
      Algorithm = 'Igg_Igg';
      r_up(1) = find(ismember(model.rxns,'ICDHyp'));
      r_up(2) = find(ismember(model.rxns,'LDH_Lm_b'));
      r_up(3) = find(ismember(model.rxns,'r1290_b'));
      r_up(4) = find(ismember(model.rxns,'FADRx_b'));
      r_up(5) = find(ismember(model.rxns,'RTOTAL3t_b'));
      r_null(1) = find(ismember(model.rxns,'PROD2'));
      r_null(2) = find(ismember(model.rxns,'MDHm_f'));
      r_null(3) = find(ismember(model.rxns,'AACOAT'));
      r_null(4) = find(ismember(model.rxns,'PCm'));
      r_null(5) = find(ismember(model.rxns,'PGCD'));
      r_null(6) = find(ismember(model.rxns,'NADHtpu_f'));
      r_null(7) = find(ismember(model.rxns,'ACACT1rm_b'));
      r_null(8) = find(ismember(model.rxns,'ALATA_L_b'));
      r_null(9) = find(ismember(model.rxns,'SUCOASm_b'));
      r_null(10) = find(ismember(model.rxns,'GLYtm_b'));
      r_null(11) = find(ismember(model.rxns,'r0422_b'));
  elseif (x == 45)
      Regulations = 225;
      Feed = 'D4';
      up = 0.9;
      Algorithm = 'Igg_Igg';
      r_up(1) = find(ismember(model.rxns,'MDHm_f'));
      r_up(2) = find(ismember(model.rxns,'PPA'));
      r_up(3) = find(ismember(model.rxns,'ACACT1r_b'));
      r_null(1) = find(ismember(model.rxns,'PROD2'));
      r_null(2) = find(ismember(model.rxns,'ACOAD8m'));
      r_null(3) = find(ismember(model.rxns,'CSNATm'));
      r_null(4) = find(ismember(model.rxns,'ACACt2m_b'));
      r_null(5) = find(ismember(model.rxns,'MMMm_b'));
      r_null(6) = find(ismember(model.rxns,'ARTPLM3_b'));
  elseif (x == 46)
      Regulations = 150;
      Feed = 'D4';
      up = 1;
      Algorithm = 'Igg_Igg';
      r_up(1) = find(ismember(model.rxns,'P5CRx'));
      r_up(2) = find(ismember(model.rxns,'PGK_f'));
      r_up(3) = find(ismember(model.rxns,'DHAPtc'));
      r_up(4) = find(ismember(model.rxns,'HMGCOAtm_f'));
      r_up(5) = find(ismember(model.rxns,'NADHtru'));
      r_up(6) = find(ismember(model.rxns,'r1109_f'));
      r_up(7) = find(ismember(model.rxns,'COAtp_b'));
      r_up(8) = find(ismember(model.rxns,'PROtm_b'));
      r_null(1) = find(ismember(model.rxns,'ASNNm'));
      r_null(2) = find(ismember(model.rxns,'MDH_f'));
      r_null(3) = find(ismember(model.rxns,'PSERT'));
      r_null(4) = find(ismember(model.rxns,'Htx_f'));
      r_null(5) = find(ismember(model.rxns,'SUCOASm_b'));
      r_null(6) = find(ismember(model.rxns,'ARTPLM2_b'));
  elseif (x == 47)
      Regulations = 200;
      Feed = 'D1';
      up = 1;
      Algorithm = 'Igg_Igg';
      r_up(1) = find(ismember(model.rxns,'RPI_f'));
      r_up(2) = find(ismember(model.rxns,'r0249_f'));
  elseif (x == 48)
      Regulations = 125;
      Feed = 'D4';
      up = 0.50;
      Algorithm = 'Bio_Igg';
      r_up(1) = find(ismember(model.rxns,'THD1m'));
      r_up(2) = find(ismember(model.rxns,'CHOLK'));
      r_up(3) = find(ismember(model.rxns,'DTTPDP'));
      r_up(4) = find(ismember(model.rxns,'ACS'));
      r_up(5) = find(ismember(model.rxns,'ADPCOAPTE'));
      r_up(6) = find(ismember(model.rxns,'ARTFR13_cho'));
      r_null(1) = find(ismember(model.rxns,'PROD2'));
      r_null(2) = find(ismember(model.rxns,'SMOX'));
      r_null(3) = find(ismember(model.rxns,'HMGCOAtm_f'));
      r_null(4) = find(ismember(model.rxns,'NADPHtxu'));
  elseif (x == 49)
      Regulations = 150;
      Feed = 'D4';
      up = 0.50;
      Algorithm = 'Bio_Igg';
      r_up(1) = find(ismember(model.rxns,'RE3301C'));
      r_up(2) = find(ismember(model.rxns,'CSm'));
      r_up(3) = find(ismember(model.rxns,'AACOAT'));
      r_up(4) = find(ismember(model.rxns,'C6DCc'));
      r_up(5) = find(ismember(model.rxns,'RTOTAL3FATPc'));
      r_null(1) = find(ismember(model.rxns,'HACD9m'));
      r_null(2) = find(ismember(model.rxns,'LEUTA'));
  elseif (x == 50)
      Regulations = 225;
      Feed = 'D4';
      up = 0.50;
      Algorithm = 'Bio_Igg';
      r_up(1) = find(ismember(model.rxns,'HMGLx'));
      r_up(2) = find(ismember(model.rxns,'PRPNCOAHYDm'));
      r_up(3) = find(ismember(model.rxns,'r0062'));
      r_up(4) = find(ismember(model.rxns,'HMGCOARc_b'));
      r_up(5) = find(ismember(model.rxns,'HMGCOAtx_b'));
      r_up(6) = find(ismember(model.rxns,'NADtx_b'));
      r_null(1) = find(ismember(model.rxns,'PROD2'));
      r_null(2) = find(ismember(model.rxns,'CSNATr_f'));
      r_null(3) = find(ismember(model.rxns,'PSERT'));
      r_null(4) = find(ismember(model.rxns,'HMGCOAtm_f'));
      r_null(5) = find(ismember(model.rxns,'NADtx_f'));
      r_null(6) = find(ismember(model.rxns,'r1290_f'));
  elseif (x == 51)
      Regulations = 200;
      Feed = 'D4';
      up = 0.60;
      Algorithm = 'Bio_Igg';
      r_up(1) = find(ismember(model.rxns,'TRDR'));
      r_up(2) = find(ismember(model.rxns,'O2tp_b'));
      r_up(3) = find(ismember(model.rxns,'RTOTAL3t_b'));
      r_null(1) = find(ismember(model.rxns,'PIt2m_f'));
      r_null(2) = find(ismember(model.rxns,'34HPPte_r'));
      r_null(3) = find(ismember(model.rxns,'NADPtxu'));
      r_null(4) = find(ismember(model.rxns,'O2tp_f'));
      r_null(5) = find(ismember(model.rxns,'r0202_b'));
  elseif (x == 52)
      Regulations = 225;
      Feed = 'D4';
      up = 0.60;
      Algorithm = 'Bio_Igg';
      r_up(1) = find(ismember(model.rxns,'RTOTAL2t_b'));
      r_up(2) = find(ismember(model.rxns,'RTOTAL3t_b'));
      r_null(1) = find(ismember(model.rxns,'THRA'));
      r_null(2) = find(ismember(model.rxns,'AACOAT'));
      r_null(3) = find(ismember(model.rxns,'LEUTA'));
      r_null(4) = find(ismember(model.rxns,'PSP_L'));
      r_null(5) = find(ismember(model.rxns,'r1109_b'));
  elseif (x == 53)
      Regulations = 150;
      Feed = 'D4';
      up = 0.70;
      Algorithm = 'Bio_Igg';
      r_up(1) = find(ismember(model.rxns,'CLHCOtex'));
      r_up(2) = find(ismember(model.rxns,'NKCCt_f'));
      r_up(3) = find(ismember(model.rxns,'PIt9_f'));
      r_up(4) = find(ismember(model.rxns,'ME2'));
      r_up(5) = find(ismember(model.rxns,'Htr_f'));
      r_up(6) = find(ismember(model.rxns,'PItn_f'));
      r_up(7) = find(ismember(model.rxns,'PPItm_f'));
      r_up(8) = find(ismember(model.rxns,'r1109_f'));
      r_up(9) = find(ismember(model.rxns,'r1503'));
      r_up(10) = find(ismember(model.rxns,'KHte_b'));
      r_up(11) = find(ismember(model.rxns,'PROtm_b'));
      r_null(1) = find(ismember(model.rxns,'NAHCO3_HCLt2_f'));
      r_null(2) = find(ismember(model.rxns,'PROD2'));
      r_null(3) = find(ismember(model.rxns,'CSm'));
      r_null(4) = find(ismember(model.rxns,'RE0383C_f'));
      r_null(5) = find(ismember(model.rxns,'LDH_Lm_f'));
      r_null(6) = find(ismember(model.rxns,'KHte_f'));
      r_null(7) = find(ismember(model.rxns,'NADPHtxu'));
      r_null(8) = find(ismember(model.rxns,'PIt8_b'));
      r_null(9) = find(ismember(model.rxns,'OCOAT1m_b'));
      r_null(10) = find(ismember(model.rxns,'GHMT2r_b'));
  elseif (x == 54)
      Regulations = 150;
      Feed = 'D4';
      up = 0.80;
      Algorithm = 'Bio_Igg';
      r_up(1) = find(ismember(model.rxns,'GLU5SAtmc_f'));
      r_null(1) = find(ismember(model.rxns,'GLU5SAtmc_b'));
  elseif (x == 55)
      Regulations = 200;
      Feed = 'D4';
      up = 0.80;
      Algorithm = 'Bio_Igg';
      r_up(1) = find(ismember(model.rxns,'r0009'));
      r_up(2) = find(ismember(model.rxns,'r2410'));
      r_up(3) = find(ismember(model.rxns,'H2Otly_b'));
      r_null(1) = find(ismember(model.rxns,'GLXO1'));
      r_null(2) = find(ismember(model.rxns,'PROD2'));
      r_null(3) = find(ismember(model.rxns,'ACACT1r_b'));
  elseif (x == 56)
      Regulations = 125;
      Feed = 'D4';
      up = 0.90;
      Algorithm = 'Bio_Igg';
      r_up(1) = find(ismember(model.rxns,'r1291_f'));
      r_up(2) = find(ismember(model.rxns,'r1503'));
      r_up(3) = find(ismember(model.rxns,'ACt2m_b'));
      r_null(1) = find(ismember(model.rxns,'PROD2'));
      r_null(2) = find(ismember(model.rxns,'COAtm_b'));
  elseif (x == 57)
      Regulations = 150;
      Feed = 'D4';
      up = 0.90;
      Algorithm = 'Bio_Igg';
      r_up(1) = find(ismember(model.rxns,'DNDPt40m'));
      r_up(2) = find(ismember(model.rxns,'PPAm'));
      r_up(3) = find(ismember(model.rxns,'C181CPT1'));
      r_up(4) = find(ismember(model.rxns,'COAtm_b'));
      r_up(5) = find(ismember(model.rxns,'RTOTAL3t_b'));
      r_null(1) = find(ismember(model.rxns,'PROD2'));
      r_null(2) = find(ismember(model.rxns,'r1154_b'));
  elseif (x == 58)
      Regulations = 175;
      Feed = 'D4';
      up = 0.90;
      Algorithm = 'Bio_Igg';
      r_up(1) = find(ismember(model.rxns,'G5SADr_f'));
      r_null(1) = find(ismember(model.rxns,'G5SADr_b'));
  elseif (x == 59)
      Regulations = 100;
      Feed = 'D4';
      up = 1;
      Algorithm = 'Bio_Igg';
      r_up(1) = find(ismember(model.rxns,'P5CRx'));
      r_up(2) = find(ismember(model.rxns,'L_LACtm_f'));
      r_up(3) = find(ismember(model.rxns,'ACACT1rm_f'));
      r_up(4) = find(ismember(model.rxns,'GLYC3PFADm'));
      r_up(5) = find(ismember(model.rxns,'H2Otp_f'));
      r_up(6) = find(ismember(model.rxns,'PPItm_f'));
      r_up(7) = find(ismember(model.rxns,'r0321'));
      r_up(8) = find(ismember(model.rxns,'r1291_f'));
      r_up(9) = find(ismember(model.rxns,'MDH_b'));
      r_null(1) = find(ismember(model.rxns,'L_LACtm_b'));
      r_null(2) = find(ismember(model.rxns,'PPItr_b'));
      r_null(3) = find(ismember(model.rxns,'r2521_b'));
  elseif (x == 60)
      Regulations = 150;
      Feed = 'D4';
      up = 1;
      Algorithm = 'Bio_Igg';
      r_up(1) = find(ismember(model.rxns,'P5CRx'));
      r_up(2) = find(ismember(model.rxns,'SUCOASm_f'));
      r_up(3) = find(ismember(model.rxns,'Htr_f'));
      r_up(4) = find(ismember(model.rxns,'PPItx_f'));
      r_up(5) = find(ismember(model.rxns,'NRVNCCOAtx_b'));
      r_null(1) = find(ismember(model.rxns,'GAPD_f'));
      r_null(2) = find(ismember(model.rxns,'peplys_synthesis'));
      r_null(3) = find(ismember(model.rxns,'PYRt2m_b'));
  elseif (x == 61)
      Regulations = 200;
      Feed = 'D4';
      up = 1;
      Algorithm = 'Bio_Igg';
      r_up(1) = find(ismember(model.rxns,'ACACt2m_f'));
      r_up(2) = find(ismember(model.rxns,'P5CRm'));
      r_up(3) = find(ismember(model.rxns,'PYRt2p_f'));
      r_up(4) = find(ismember(model.rxns,'GLYC3PFADm'));
      r_up(5) = find(ismember(model.rxns,'MMTSADm'));
      r_up(6) = find(ismember(model.rxns,'r1464_f'));
      r_up(7) = find(ismember(model.rxns,'INSTt2r_b'));
      r_up(8) = find(ismember(model.rxns,'ACACT1r_b'));
      r_up(9) = find(ismember(model.rxns,'Htx_b'));
      r_null(1) = find(ismember(model.rxns,'INSTt2r_f'));
      r_null(2) = find(ismember(model.rxns,'GAPD_f'));
      r_null(3) = find(ismember(model.rxns,'PPAm'));
      r_null(4) = find(ismember(model.rxns,'r0009'));
      r_null(5) = find(ismember(model.rxns,'MDH_f'));
      r_null(6) = find(ismember(model.rxns,'PGCD'));
      r_null(7) = find(ismember(model.rxns,'CO2tp_f'));
      r_null(8) = find(ismember(model.rxns,'FORt2m_r'));
      r_null(9) = find(ismember(model.rxns,'LACLt'));
      r_null(10) = find(ismember(model.rxns,'RPE_b'));
  elseif (x == 62)
      Regulations = 200;
      Feed = 'D4';
      up = 1;
      Algorithm = 'Bio_Igg';
      r_up(1) = find(ismember(model.rxns,'P5CR'));
      r_up(2) = find(ismember(model.rxns,'CYTDK1'));
      r_up(3) = find(ismember(model.rxns,'DNDPt47m'));
      r_up(4) = find(ismember(model.rxns,'ACACT1rm_f'));
      r_up(5) = find(ismember(model.rxns,'GLYC3PFADm'));
      r_up(6) = find(ismember(model.rxns,'r0556_r'));
      r_up(7) = find(ismember(model.rxns,'COAtm_b'));
      r_up(8) = find(ismember(model.rxns,'r1290_b'));
      r_null(1) = find(ismember(model.rxns,'SOAT11_cho'));
      r_null(2) = find(ismember(model.rxns,'NADPtxu'));
      r_null(3) = find(ismember(model.rxns,'NADtx_f'));
      r_null(4) = find(ismember(model.rxns,'PPItm_f'));
      r_null(5) = find(ismember(model.rxns,'r0173_f'));
      r_null(6) = find(ismember(model.rxns,'PHETHPTOX'));
      r_null(7) = find(ismember(model.rxns,'r0388_b'));
  elseif (x == 63)
      Regulations = 225;
      Feed = 'D4';
      up = 1;
      Algorithm = 'Bio_Igg';
      r_up(1) = find(ismember(model.rxns,'G3PD2m'));
      r_up(2) = find(ismember(model.rxns,'Htr_f'));
      r_up(3) = find(ismember(model.rxns,'MMTSADm'));
      r_up(4) = find(ismember(model.rxns,'r1109_f'));
      r_up(5) = find(ismember(model.rxns,'r1291_f'));
      r_up(6) = find(ismember(model.rxns,'ALATA_L_b'));
      r_up(7) = find(ismember(model.rxns,'COAtp_b'));
      r_up(8) = find(ismember(model.rxns,'NADHtpu_b'));
      r_null(1) = find(ismember(model.rxns,'PROD2'));
      r_null(2) = find(ismember(model.rxns,'PIt2m_f'));
      r_null(3) = find(ismember(model.rxns,'PYRtm_f'));
      r_null(4) = find(ismember(model.rxns,'HMGCOASi'));
      r_null(5) = find(ismember(model.rxns,'HMGCOARx'));
      r_null(6) = find(ismember(model.rxns,'r1434'));
      r_null(7) = find(ismember(model.rxns,'HMGCOAtm_b'));
      r_null(8) = find(ismember(model.rxns,'Htx_b'));
  elseif (x == 64)
      Regulations = 50;
      Feed = 'D4';
      up = 1.1;
      Algorithm = 'Bio_Igg';
      r_null(1) = find(ismember(model.rxns,'RNDR1'));
      r_null(2) = find(ismember(model.rxns,'PIt2m_f'));
      r_null(3) = find(ismember(model.rxns,'GAPD_f'));
      r_null(4) = find(ismember(model.rxns,'ACACT1rm_f'));
      r_null(5) = find(ismember(model.rxns,'ASNNm'));
      r_null(6) = find(ismember(model.rxns,'FUMm_f'));
      r_null(7) = find(ismember(model.rxns,'PEPCKm'));
      r_null(8) = find(ismember(model.rxns,'r1148_f'));
      r_null(9) = find(ismember(model.rxns,'SPMTDe_r'));
      r_null(10) = find(ismember(model.rxns,'peplys_synthesis'));
      r_null(11) = find(ismember(model.rxns,'FADRx_b'));
  elseif (x == 65)
      Regulations = 175;
      Feed = 'D3';
      up = 0.50;
      Algorithm = 'Igg_Igg';
      r_up(1) = find(ismember(model.rxns,'SO4OXAtex2'));
      r_up(2) = find(ismember(model.rxns,'ACTNMO'));
      r_up(3) = find(ismember(model.rxns,'r0642_r'));
      r_up(4) = find(ismember(model.rxns,'r1154_b'));
      r_null(1) = find(ismember(model.rxns,'VALTA'));
      r_null(2) = find(ismember(model.rxns,'r0084_b'));
  elseif (x == 66)
      Regulations = 100;
      Feed = 'D3';
      up = 0.90;
      Algorithm = 'Igg_Igg';
      r_up(1) = find(ismember(model.rxns,'GLUt2m_f'));
      r_null(1) = find(ismember(model.rxns,'PROD2'));
      r_null(2) = find(ismember(model.rxns,'PIt2m_f'));
      r_null(3) = find(ismember(model.rxns,'COAtm_f'));
      r_null(4) = find(ismember(model.rxns,'PPItm_f'));
  elseif (x == 67)
      Regulations = 175;
      Feed = 'D3';
      up = 0.90;
      Algorithm = 'Igg_Igg';
      r_up(1) = find(ismember(model.rxns,'NMNtn'));
      r_up(2) = find(ismember(model.rxns,'PIK4n_cho'));
      r_up(3) = find(ismember(model.rxns,'r0383'));
      r_up(4) = find(ismember(model.rxns,'r1290_f'));
      r_up(5) = find(ismember(model.rxns,'CSNATp_b'));
      r_null(1) = find(ismember(model.rxns,'GLU5SAtmc_b'));
      r_null(2) = find(ismember(model.rxns,'r1290_b'));
  elseif (x == 68)
      Regulations = 225;
      Feed = 'D3';
      up = 0.90;
      Algorithm = 'Igg_Igg';
      r_up(1) = find(ismember(model.rxns,'r1454'));
      r_up(2) = find(ismember(model.rxns,'r1459_b'));
      r_null(1) = find(ismember(model.rxns,'ACt2m_f'));
      r_null(2) = find(ismember(model.rxns,'CSm'));
      r_null(3) = find(ismember(model.rxns,'PPM_f'));
      r_null(4) = find(ismember(model.rxns,'ALDD2xm'));
      r_null(5) = find(ismember(model.rxns,'G6PDH2rer_f'));
      r_null(6) = find(ismember(model.rxns,'r0062'));
      r_null(7) = find(ismember(model.rxns,'r1290_f'));
      r_null(8) = find(ismember(model.rxns,'G5SADr_b'));
      r_null(9) = find(ismember(model.rxns,'MDHm_b'));
      r_null(10) = find(ismember(model.rxns,'PPItx_b'));
  elseif (x == 69)
      Regulations = 175;
      Feed = 'D3';
      up = 1.1;
      Algorithm = 'Igg_Igg';
      r_null(1) = find(ismember(model.rxns,'SPHK21c'));
      r_null(2) = find(ismember(model.rxns,'ASNN'));
      r_null(3) = find(ismember(model.rxns,'CEPTE_cho'));
      r_null(4) = find(ismember(model.rxns,'r0060'));
      r_null(5) = find(ismember(model.rxns,'2AMACHYD'));
      r_null(6) = find(ismember(model.rxns,'ACACT1rm_f'));
      r_null(7) = find(ismember(model.rxns,'ENO_f'));
      r_null(8) = find(ismember(model.rxns,'MCCCrm'));
      r_null(9) = find(ismember(model.rxns,'MI34PP'));
      r_null(10) = find(ismember(model.rxns,'ASPGLUm_f'));
      r_null(11) = find(ismember(model.rxns,'PSDm_cho'));
      r_null(12) = find(ismember(model.rxns,'2AMACSULT'));
      r_null(13) = find(ismember(model.rxns,'34HPPte_r'));
      r_null(14) = find(ismember(model.rxns,'ARTPLM1_f'));
      r_null(15) = find(ismember(model.rxns,'NCAMDe_r'));
      r_null(16) = find(ismember(model.rxns,'GMPS'));
      r_null(17) = find(ismember(model.rxns,'SPRMTDe'));
      r_null(18) = find(ismember(model.rxns,'G5SADr_b'));
      r_null(19) = find(ismember(model.rxns,'ASPTA_b'));
  elseif (x == 70)
      Regulations = 100;
      Feed = 'D3';
      up = 0.50;
      Algorithm = 'Bio_Igg';
      r_up(1) = find(ismember(model.rxns,'GLYKm'));
      r_up(2) = find(ismember(model.rxns,'PGK_f'));
      r_up(3) = find(ismember(model.rxns,'FADH2ETC'));
      r_up(4) = find(ismember(model.rxns,'r0383'));
      r_up(5) = find(ismember(model.rxns,'FUMm_b'));
      r_null(1) = find(ismember(model.rxns,'HISD'));
      r_null(2) = find(ismember(model.rxns,'MDH_f'));
      r_null(3) = find(ismember(model.rxns,'34HPPte_r'));
  elseif (x == 71)
      Regulations = 225;
      Feed = 'D2';
      up = 0.50;
      Algorithm = 'Igg_Igg';
      r_up(1) = find(ismember(model.rxns,'ACS'));
      r_up(2) = find(ismember(model.rxns,'DHAPtc'));
      r_up(3) = find(ismember(model.rxns,'r1525'));
      r_up(4) = find(ismember(model.rxns,'RE1530C_f'));
      r_up(5) = find(ismember(model.rxns,'H2Otn_b'));
      r_up(6) = find(ismember(model.rxns,'PE_CHOtm_b'));
      r_null(1) = find(ismember(model.rxns,'PAFS_cho'));
      r_null(2) = find(ismember(model.rxns,'LEUTA'));
      r_null(3) = find(ismember(model.rxns,'r1290_f'));
      r_null(4) = find(ismember(model.rxns,'RE1530C_b'));
  elseif (x == 72)
      Regulations = 225;
      Feed = 'D2';
      up = 0.50;
      Algorithm = 'Igg_Igg';
      r_up(1) = find(ismember(model.rxns,'RPE_f'));
      r_null(1) = find(ismember(model.rxns,'G6PDH2r'));
      r_null(2) = find(ismember(model.rxns,'GNDer'));
      r_null(3) = find(ismember(model.rxns,'TKT2_b'));
  elseif (x == 73)
      Regulations = 150;
      Feed = 'D2';
      up = 0.90;
      Algorithm = 'Igg_Igg';
      r_up(1) = find(ismember(model.rxns,'PIt2m_f'));
      r_up(2) = find(ismember(model.rxns,'r0086'));
      r_up(3) = find(ismember(model.rxns,'r2420_f'));
      r_up(4) = find(ismember(model.rxns,'G5SADr_b'));
      r_up(5) = find(ismember(model.rxns,'H2Otly_b'));
      r_null(1) = find(ismember(model.rxns,'LEUTA'));
      r_null(2) = find(ismember(model.rxns,'r2472'));
      r_null(3) = find(ismember(model.rxns,'ACt2m_b'));
      r_null(4) = find(ismember(model.rxns,'PYRt2m_b'));
      r_null(5) = find(ismember(model.rxns,'r1109_b'));
      r_null(6) = find(ismember(model.rxns,'r2420_b'));
  elseif (x == 74)
      Regulations = 25;
      Feed = 'D2';
      up = 1.1;
      Algorithm = 'Igg_Igg';
      r_null(1) = find(ismember(model.rxns,'PROD2'));
      r_null(2) = find(ismember(model.rxns,'PIt2m_f'));
      r_null(3) = find(ismember(model.rxns,'GHMT2rm_f'));
      r_null(4) = find(ismember(model.rxns,'LEUTA'));
      r_null(5) = find(ismember(model.rxns,'PROtm_f'));
      r_null(6) = find(ismember(model.rxns,'r1464_f'));
      r_null(7) = find(ismember(model.rxns,'L_LACtm_b'));
      r_null(8) = find(ismember(model.rxns,'ORNt4m_b'));
      r_null(9) = find(ismember(model.rxns,'ASPGLUm_b'));
      r_null(10) = find(ismember(model.rxns,'PYRt2m_b'));
  elseif (x == 75)
      Regulations = 175;
      Feed = 'D2';
      up = 1.1;
      Algorithm = 'Igg_Igg';
      r_null(1) = find(ismember(model.rxns,'FAS160ACP'));
      r_null(2) = find(ismember(model.rxns,'ARD1_cho'));
      r_null(3) = find(ismember(model.rxns,'PROD2'));
      r_null(4) = find(ismember(model.rxns,'THRA'));
      r_null(5) = find(ismember(model.rxns,'HACD9m'));
      r_null(6) = find(ismember(model.rxns,'MTHFC_f'));
      r_null(7) = find(ismember(model.rxns,'PYRtm_f'));
      r_null(8) = find(ismember(model.rxns,'CSm'));
      r_null(9) = find(ismember(model.rxns,'ACOAO7p'));
      r_null(10) = find(ismember(model.rxns,'PPAm'));
      r_null(11) = find(ismember(model.rxns,'GHMT2rm_f'));
      r_null(12) = find(ismember(model.rxns,'LEUTA'));
      r_null(13) = find(ismember(model.rxns,'NNMT'));
      r_null(14) = find(ismember(model.rxns,'5MTHFt2_r'));
      r_null(15) = find(ismember(model.rxns,'PSDm_cho'));
      r_null(16) = find(ismember(model.rxns,'3MOBt2im'));
      r_null(17) = find(ismember(model.rxns,'DOLP_ter'));
      r_null(18) = find(ismember(model.rxns,'r1290_f'));
      r_null(19) = find(ismember(model.rxns,'RE3242C'));
      r_null(20) = find(ismember(model.rxns,'OCOAT1m_b'));
      r_null(21) = find(ismember(model.rxns,'FUM_b'));
      r_null(22) = find(ismember(model.rxns,'RPE_b'));
      r_null(23) = find(ismember(model.rxns,'r0422_b'));
      r_null(24) = find(ismember(model.rxns,'r1109_b'));
  elseif (x == 76)
      Regulations = 225;
      Feed = 'D2';
      up = 0.50;
      Algorithm = 'Bio_Igg';
      r_up(1) = find(ismember(model.rxns,'r0086'));
      r_up(2) = find(ismember(model.rxns,'ACS'));
      r_up(3) = find(ismember(model.rxns,'IMPD'));
      r_up(4) = find(ismember(model.rxns,'MCDm'));
      r_up(5) = find(ismember(model.rxns,'MMSAD3m'));
      r_up(6) = find(ismember(model.rxns,'HMGCOAtm_f'));
      r_up(7) = find(ismember(model.rxns,'THR3D'));
      r_null(1) = find(ismember(model.rxns,'KAS8'));
      r_null(2) = find(ismember(model.rxns,'VALTA'));
      r_null(3) = find(ismember(model.rxns,'G5SADr_b'));
      r_null(4) = find(ismember(model.rxns,'TKT2_b'));
      r_null(5) = find(ismember(model.rxns,'HMGCOAtm_b'));
  elseif (x == 77)
      Regulations = 175;
      Feed = 'D2';
      up = 0.60;
      Algorithm = 'Bio_Igg';
      r_up(1) = find(ismember(model.rxns,'LDH_D_f'));
      r_up(2) = find(ismember(model.rxns,'BIOCYTtn'));
      r_up(3) = find(ismember(model.rxns,'r0062'));
      r_up(4) = find(ismember(model.rxns,'r0357'));
      r_up(5) = find(ismember(model.rxns,'RTOT6_cho'));
      r_null(1) = find(ismember(model.rxns,'UGCGc'));
      r_null(2) = find(ismember(model.rxns,'PROD2'));
      r_null(3) = find(ismember(model.rxns,'OCOAT1m_f'));
      r_null(4) = find(ismember(model.rxns,'ASNN'));
      r_null(5) = find(ismember(model.rxns,'PPAm'));
      r_null(6) = find(ismember(model.rxns,'ALR3'));
      r_null(7) = find(ismember(model.rxns,'CLOXAtex2'));
      r_null(8) = find(ismember(model.rxns,'ASNtm'));
      r_null(9) = find(ismember(model.rxns,'FADH2tru'));
      r_null(10) = find(ismember(model.rxns,'GLNtm'));
      r_null(11) = find(ismember(model.rxns,'NADHtru'));
      r_null(12) = find(ismember(model.rxns,'NADPtru'));
      r_null(13) = find(ismember(model.rxns,'r0431'));
      r_null(14) = find(ismember(model.rxns,'r1434'));
  elseif (x == 78)
      Regulations = 225;
      Feed = 'D1';
      up = 0.50;
      Algorithm = 'Bio_Igg';
      r_up(1) = find(ismember(model.rxns,'PI34P3P'));
      r_up(2) = find(ismember(model.rxns,'DNDPt60m'));
      r_up(3) = find(ismember(model.rxns,'DRBK'));
      r_up(4) = find(ismember(model.rxns,'GND'));
      r_up(5) = find(ismember(model.rxns,'PYDXtr_r'));
      r_up(6) = find(ismember(model.rxns,'r0221'));
      r_up(7) = find(ismember(model.rxns,'r0366'));
      r_null(1) = find(ismember(model.rxns,'GLXO1'));
      r_null(2) = find(ismember(model.rxns,'ABUTt2r_r'));
      r_null(3) = find(ismember(model.rxns,'34HPPte_r'));
  elseif (x == 79)
      Regulations = 225;
      Feed = 'D1';
      up = 0.50;
      Algorithm = 'Bio_Igg';
      r_up(1) = find(ismember(model.rxns,'BILGLCURte'));
      r_up(2) = find(ismember(model.rxns,'THD1m'));
      r_up(3) = find(ismember(model.rxns,'MMCDm'));
      r_up(4) = find(ismember(model.rxns,'r0086'));
      r_up(5) = find(ismember(model.rxns,'RTOTAL3t_b'));
      r_null(1) = find(ismember(model.rxns,'MI4PK'));
      r_null(2) = find(ismember(model.rxns,'PPAer'));
      r_null(3) = find(ismember(model.rxns,'PPAm'));
      r_null(4) = find(ismember(model.rxns,'ACETONEt2m_r'));
      r_null(5) = find(ismember(model.rxns,'PGCD'));
      r_null(6) = find(ismember(model.rxns,'34HPPte_r'));
  elseif (x == 80)
      Regulations = 225;
      Feed = 'D1';
      up = 0.50;
      Algorithm = 'Bio_Igg';
      r_up(1) = find(ismember(model.rxns,'RPE_f'));
      r_null(1) = find(ismember(model.rxns,'GND'));
      r_null(2) = find(ismember(model.rxns,'GNDer'));
      r_null(3) = find(ismember(model.rxns,'TKT2_b'));
  elseif (x == 81)
      Regulations = 100;
      Feed = 'D1';
      up = 0.60;
      Algorithm = 'Bio_Igg';
      r_up(1) = find(ismember(model.rxns,'ACACT1rm_f'));
      r_up(2) = find(ismember(model.rxns,'ACS'));
      r_up(3) = find(ismember(model.rxns,'PYRt2m_f'));
      r_up(4) = find(ismember(model.rxns,'PItx_r'));
      r_up(5) = find(ismember(model.rxns,'PPItm_f'));
      r_null(1) = find(ismember(model.rxns,'GAPD_f'));
      r_null(2) = find(ismember(model.rxns,'34HPPte_r'));
      r_null(3) = find(ismember(model.rxns,'NADPtxu'));
      r_null(4) = find(ismember(model.rxns,'r1290_f'));
  elseif (x == 82)
      Regulations = 175;
      Feed = 'D1';
      up = 0.60;
      Algorithm = 'Bio_Igg';
      r_up(1) = find(ismember(model.rxns,'GLNS'));
      r_up(2) = find(ismember(model.rxns,'NMNS'));
      r_up(3) = find(ismember(model.rxns,'DNDPt39m'));
      r_up(4) = find(ismember(model.rxns,'PYDXNK'));
      r_up(5) = find(ismember(model.rxns,'BTND1'));
      r_up(6) = find(ismember(model.rxns,'GK1'));
      r_up(7) = find(ismember(model.rxns,'r1528'));
      r_up(8) = find(ismember(model.rxns,'PEt_cho_r'));
      r_up(9) = find(ismember(model.rxns,'PIK4n_cho'));
      r_up(10) = find(ismember(model.rxns,'PRISTCOAtcx'));
      r_up(11) = find(ismember(model.rxns,'VITKtl'));
      r_null(1) = find(ismember(model.rxns,'CHOLtu_f'));
      r_null(2) = find(ismember(model.rxns,'ASNNm'));
      r_null(3) = find(ismember(model.rxns,'34HPPte_r'));
      r_null(4) = find(ismember(model.rxns,'r1148_f'));
  elseif (x == 83)
      Regulations = 100;
      Feed = 'D1';
      up = 0.70;
      Algorithm = 'Bio_Igg';
      r_up(1) = find(ismember(model.rxns,'PGK_f'));
      r_up(2) = find(ismember(model.rxns,'GALU'));
      r_up(3) = find(ismember(model.rxns,'P5CDm'));
      r_up(4) = find(ismember(model.rxns,'r0062'));
      r_up(5) = find(ismember(model.rxns,'r1144'));
      r_null(1) = find(ismember(model.rxns,'CHOLtu_f'));
      r_null(2) = find(ismember(model.rxns,'34HPPte_r'));
      r_null(3) = find(ismember(model.rxns,'FADRx_b'));
  elseif (x == 84)
      Regulations = 75;
      Feed = 'D1';
      up = 0.90;
      Algorithm = 'Bio_Igg';
      r_up(1) = find(ismember(model.rxns,'PIt2m_f'));
      r_up(2) = find(ismember(model.rxns,'O2tp_f'));
      r_up(3) = find(ismember(model.rxns,'r1501'));
      r_up(4) = find(ismember(model.rxns,'H2Otly_b'));
      r_null(1) = find(ismember(model.rxns,'AASAD3_r'));
      r_null(2) = find(ismember(model.rxns,'MDH_f'));
      r_null(3) = find(ismember(model.rxns,'34HPPte_r'));
      r_null(4) = find(ismember(model.rxns,'HMGCOAtm_b'));
  elseif (x == 85)
      Regulations = 200;
      Feed = 'D1';
      up = 1;
      Algorithm = 'Bio_Igg';
      r_up(1) = find(ismember(model.rxns,'RPI_f'));
      r_up(2) = find(ismember(model.rxns,'r0249_f'));
      r_null(1) = find(ismember(model.rxns,'TKT2_b'));
  elseif (x == 86)
      Regulations = 25;
      Feed = 'D1';
      up = 1.1;
      Algorithm = 'Bio_Igg';
      r_null(1) = find(ismember(model.rxns,'GLXO1'));
      r_null(2) = find(ismember(model.rxns,'PIt2m_f'));
      r_null(3) = find(ismember(model.rxns,'34HPPte_r'));
      r_null(4) = find(ismember(model.rxns,'GLUt2m_b'));
      r_null(5) = find(ismember(model.rxns,'SUCOASm_b'));
  elseif (x == 87)
      Regulations = 150;
      Feed = 'D1';
      up = 1.1;
      Algorithm = 'Bio_Igg';
      r_null(1) = find(ismember(model.rxns,'PYRtm_f'));
      r_null(2) = find(ismember(model.rxns,'ACACT1rm_f'));
      r_null(3) = find(ismember(model.rxns,'PUNP3_f'));
      r_null(4) = find(ismember(model.rxns,'34HPPte_r'));
      r_null(5) = find(ismember(model.rxns,'LYStn_f'));
      r_null(6) = find(ismember(model.rxns,'GLUt2m_b'));
      r_null(7) = find(ismember(model.rxns,'H2Otn_b'));
      r_null(8) = find(ismember(model.rxns,'PPItx_b'));
      r_null(9) = find(ismember(model.rxns,'r0422_b'));
      r_null(10) = find(ismember(model.rxns,'r0424_b'));
  elseif (x == 88)
      Regulations = 150;
      Feed = 'D1';
      up = 1.1;
      Algorithm = 'Bio_Igg';
      r_null(1) = find(ismember(model.rxns,'THRA'));
      r_null(2) = find(ismember(model.rxns,'PPAm'));
      r_null(3) = find(ismember(model.rxns,'AACOAT'));
      r_null(4) = find(ismember(model.rxns,'ALR3'));
      r_null(5) = find(ismember(model.rxns,'HMGCOASim'));
      r_null(6) = find(ismember(model.rxns,'MI14PP'));
      r_null(7) = find(ismember(model.rxns,'PCm'));
      r_null(8) = find(ismember(model.rxns,'34HPPte_r'));
      r_null(9) = find(ismember(model.rxns,'ARTPLM1_f'));
      r_null(10) = find(ismember(model.rxns,'FADH2tru'));
      r_null(11) = find(ismember(model.rxns,'Q10H2e'));
      r_null(12) = find(ismember(model.rxns,'SPRMTDe'));
      r_null(13) = find(ismember(model.rxns,'ASPTA_b'));
      r_null(14) = find(ismember(model.rxns,'GLUt2m_b'));
      r_null(15) = find(ismember(model.rxns,'L_LACtm_b'));
      r_null(16) = find(ismember(model.rxns,'MDH_b'));
      r_null(17) = find(ismember(model.rxns,'SUCOASm_b'));
      r_null(18) = find(ismember(model.rxns,'PYRt2m_b'));
      r_null(19) = find(ismember(model.rxns,'r0422_b'));
      r_null(20) = find(ismember(model.rxns,'FADRx_b'));
  elseif (x == 89)
      Regulations = 100;
      Feed = 'D2';
      up = 0.4;
      Algorithm = 'Bio_Igg';
      r_up(1) = find(ismember(model.rxns,'GTHDH'));
      r_up(2) = find(ismember(model.rxns,'THFt2_r'));
      r_up(3) = find(ismember(model.rxns,'R_group_phosphotase_1'));
      r_up(4) = find(ismember(model.rxns,'TAGt_cho_r'));
      r_up(5) = find(ismember(model.rxns,'RTOTAL2t_b'));
      r_null(1) = find(ismember(model.rxns,'ASNN'));
      r_null(2) = find(ismember(model.rxns,'PPAer'));
      r_null(3) = find(ismember(model.rxns,'ASNNm'));
      r_null(4) = find(ismember(model.rxns,'COAtm_f'));
      r_null(5) = find(ismember(model.rxns,'ARTPLM1_f'));
      r_null(6) = find(ismember(model.rxns,'r1291_f'));
      r_null(7) = find(ismember(model.rxns,'r1454'));
      r_null(8) = find(ismember(model.rxns,'RTOTAL2t_f'));




    end
else
    if (x == 1)
        Algorithm = 'Igg_Igg';
        Feed = 'D2';
        up = 0.8;
        r_up(1) = find(ismember(model.rxns,'GLCNACT_cho'));
        r_null(1) = find(ismember(model.rxns,'ACOAD10m'));
        r_null(2) = find(ismember(model.rxns,'ACGAMK'));
        r_null(3) = find(ismember(model.rxns,'3MOPte_r'));
        r_null(4) = find(ismember(model.rxns,'GALtly'));
        r_null(5) = find(ismember(model.rxns,'FBA_b'));
        
        % 1. 'GLCNACT_cho': 'naglc2p[c] + uacgam[c]  -> h[c] + chito2pdol[c] + udp[c] '
        % 1. 'ACOAD10m': 'fad[m] + 2mbcoa[m]  -> fadh2[m] + 2mb2coa[m] '
        % 2. 'ACGAMK': 'atp[c] + acgam[c]  -> h[c] + adp[c] + acgam6p[c] '
        % 3. '3MOPte': '3mop[e]  <=> 3mop[c] '
        % 4. 'ACGALtly': 'acgal[l]  -> acgal[c] '
        % 5 .'FBA': 'fdp[c]  <=> dhap[c] + dhap[c] 
        
    elseif (x == 2)
        Algorithm = 'Igg_Igg';
        Feed = 'D4';
        up = 0.5;       
        r_null(1) = find(ismember(model.rxns,'THBPT4ACAMDASE'));
        r_null(2) = find(ismember(model.rxns,'TRPO2'));
        r_null(3) = find(ismember(model.rxns,'CBPSam'));       
        
    elseif (x == 3)
        Algorithm = 'Igg_Igg';
        Feed = 'D4';
        up = 0.9;
        r_null(1) = find(ismember(model.rxns,'KYN3OX'));
        r_null(2) = find(ismember(model.rxns,'r0239'));
        r_null(3) = find(ismember(model.rxns,'KYN'));
        r_null(4) = find(ismember(model.rxns,'TRPHYDRO2'));

	    % 'KYN3OX': 'o2[c] + h[c] + nadph[c] + Lkynr[c]  -> nadp[c] + h2o[c] + hLkynr[c] '
        % 'TRPHYDRO2': 'o2[c] + trp_L[c] + thbpt[c]  -> 5htrp[c] + thbpt4acam[c] '
        % 'ANTHte': 'anth[e]  <=> anth[c] '

        
    elseif (x == 4)
        Algorithm = 'Bio_Igg';
        Feed = 'D4';
        up = 1.1;
        r_null(1) = find(ismember(model.rxns,'THBPT4ACAMDASE'));
        r_null(2) = find(ismember(model.rxns,'TRPO2'));
        r_null(3) = find(ismember(model.rxns,'UGALNACtg'));    
        % 'THBPT4ACAMDASE': 'thbpt4acam[c]  -> h2o[c] + dhbpt[c] '
        % 'TRPO2': 'o2[c] + trp_L[c]  -> Lfmkynr[c] '
        % 'UGALNACtg': 'ump[g] + udpacgal[c]  -> udpacgal[g] + ump[c] 
        
    elseif (x == 5)
        Algorithm = 'Igg_Igg';
        Feed = 'D2';
        up = 0.8;
        r_up(1) = find(ismember(model.rxns,'MANtg'));      
        r_null(1) = find(ismember(model.rxns,'HIBDm'));
        r_null(2) = find(ismember(model.rxns,'3MOBte_r'));
        r_null(3) = find(ismember(model.rxns,'ACGAMtly'));
        r_null(4) = find(ismember(model.rxns,'FBA_b'));
        r_null(5) = find(ismember(model.rxns,'TALA_b'));
        
    elseif (x == 6)
        Algorithm = 'Bio_Igg';
        Feed = 'D1';
        up = 1.1;
        r_null(1) = find(ismember(model.rxns,'PRDXl'));
        r_null(2) = find(ismember(model.rxns,'TRPO2'));
        r_null(3) = find(ismember(model.rxns,'TRPHYDRO2'));
        r_null(4) = find(ismember(model.rxns,'UAG2EMAi'));
        r_null(5) = find(ismember(model.rxns,'H2Otly_f'));
        
        % 1. 'PRDXl': 'h2o2[l] + meoh[l]  -> 2 h2o[l] + fald[l] '
        % 2. 'TRPO2': 'o2[c] + trp_L[c]  -> Lfmkynr[c] '
        % 3. 'TRPHYDRO2': 'o2[c] + trp_L[c] + thbpt[c]  -> 5htrp[c] + thbpt4acam[c] '
        % 4 .'UAG2EMAi': ''h2o[c] + uacgam[c]  -> h[c] + udp[c] + acmana[c] '
        % 5. 'H2Otly': ''h2o[c]  <=> h2o[l] '
       
                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% NO EPO IN THE MODEL %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    elseif (x == 7)
        Algorithm = 'Bio_Igg';
        Feed = 'D4';
        up = 1.1;
        r_null(1) = find(ismember(model.rxns,'TRPO2'));
        r_null(2) = find(ismember(model.rxns,'TRPHYDRO2'));

    elseif (x == 8)
        Algorithm = 'Bio_Igg';
        Feed = 'D1';
        up = 1.1;
        r_null(1) = find(ismember(model.rxns,'AHC'));
        r_null(2) = find(ismember(model.rxns,'CEPTE_cho'));
        r_null(3) = find(ismember(model.rxns,'CHOLK'));
        r_null(4) = find(ismember(model.rxns,'MGCHrm'));
        r_null(5) = find(ismember(model.rxns,'4MOPte_r'));


    end
end

if ~exist('r_up','var')
    r_up = [];
end
if ~exist('r_null','var')
    r_null = [];
end

end