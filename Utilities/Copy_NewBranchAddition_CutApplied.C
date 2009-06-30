void copy() {

   TChain *oldtree = new TChain("analysis") ;
   oldtree->Add("TauolaTTbar_S08ReDigi.root");
  
   cout<<oldtree->GetMaxTreeSize()<<endl;
   long long max_tree_size1 = 200000000000LL; // 200 GB
   if(oldtree->GetMaxTreeSize() < max_tree_size1 ) {
     oldtree->SetMaxTreeSize(max_tree_size1);
   }

   cout<<oldtree->GetMaxTreeSize()<<endl;

   double  nEvents  = (int) oldtree->GetEntries();
   double Xsection = 241.7;
   double FilterEff= 1.0;
   double scale = Xsection*100*FilterEff/nEvents;
   cout<<nEvents<<"        "<<Xsection<<"     "<<FilterEff<<"   "<<scale<<endl;
   // -- Disable all branches ...
   oldtree->SetBranchStatus("*",0);
   

   // -- ... and switch on those you'd like to write out into the new tree


   oldtree->SetBranchStatus("nt_overflow",1);// &nt_overflow, &b_nt_overflow);
   oldtree->SetBranchStatus("evt_weight",1);// &evt_weight, &b_evt_weight);
   oldtree->SetBranchStatus("proc_id",1);// &proc_id, &b_proc_id);
   oldtree->SetBranchStatus("genEventScale",1);// &genEventScale, &b_genEventScale);
   oldtree->SetBranchStatus("filter_eff",1);// &filter_eff, &b_filter_eff);
   oldtree->SetBranchStatus("cross_section",1);// &cross_section, &b_cross_section);
   oldtree->SetBranchStatus("missEn_Et_pat",1);// &missEn_Et_pat, &b_missEn_Et_pat);
   oldtree->SetBranchStatus("missEn_Ex_pat",1);// &missEn_Ex_pat, &b_missEn_Ex_pat);
   oldtree->SetBranchStatus("missEn_Ey_pat",1);// &missEn_Ey_pat, &b_missEn_Ey_pat);
   oldtree->SetBranchStatus("missEn_unCorrEt_pat",1);// &missEn_unCorrEt_pat, &b_missEn_unCorrEt_pat);
   oldtree->SetBranchStatus("missEn_allCorrEx_pat",1);// &missEn_allCorrEx_pat, &b_missEn_allCorrEx_pat);
   oldtree->SetBranchStatus("missEn_allCorrEy_pat",1);// &missEn_allCorrEy_pat, &b_missEn_allCorrEy_pat);
   oldtree->SetBranchStatus("missEn_allCorrSumEt_pat",1);// &missEn_allCorrSumEt_pat, &b_missEn_allCorrSumEt_pat);
   oldtree->SetBranchStatus("missEn_jetCorrEx_pat",1);// &missEn_jetCorrEx_pat, &b_missEn_jetCorrEx_pat);
   oldtree->SetBranchStatus("missEn_jetCorrEy_pat",1);// &missEn_jetCorrEy_pat, &b_missEn_jetCorrEy_pat);
   oldtree->SetBranchStatus("missEn_jetCorrSumEt_pat",1);// &missEn_jetCorrSumEt_pat, &b_missEn_jetCorrSumEt_pat);
   oldtree->SetBranchStatus("missEn_muCorrEx_pat",1);// &missEn_muCorrEx_pat, &b_missEn_muCorrEx_pat);
   oldtree->SetBranchStatus("missEn_muCorrEy_pat",1);// &missEn_muCorrEy_pat, &b_missEn_muCorrEy_pat);
   oldtree->SetBranchStatus("missEn_muCorrSumEt_pat",1);// &missEn_muCorrSumEt_pat, &b_missEn_muCorrSumEt_pat);
   oldtree->SetBranchStatus("sumEt_pat",1);// &sumEt_pat, &b_sumEt_pat);

   oldtree->SetBranchStatus("missEn_Et_pflow",1);// &missEn_Et_pflow, &b_missEn_Et_pflow);
   oldtree->SetBranchStatus("missEn_Ex_pflow",1);// &missEn_Ex_pflow, &b_missEn_Ex_pflow);
   oldtree->SetBranchStatus("missEn_Ey_pflow",1);// &missEn_Ey_pflow, &b_missEn_Ey_pflow);
   oldtree->SetBranchStatus("sumEt_pflow",1);// &sumEt_pflow, &b_sumEt_pflow);

   oldtree->SetBranchStatus("nE_pat",1);// &nE_pat, &b_nE_pat);
   //oldtree->SetBranchStatus("e_TMbf_pat",1);// e_TMbf_pat, &b_e_TMbf_pat);
   //oldtree->SetBranchStatus("e_TM_pat",1);// e_TM_pat, &b_e_TM_pat);
  // oldtree->SetBranchStatus("e_RM_pat",1);// e_RM_pat, &b_e_RM_pat);
   oldtree->SetBranchStatus("e_P_pat",1);// e_P_pat, &b_e_P_pat);
   oldtree->SetBranchStatus("e_Pt_pat",1);// e_Pt_pat, &b_e_Pt_pat);
   oldtree->SetBranchStatus("e_Px_pat",1);// e_Px_pat, &b_e_Px_pat);
   oldtree->SetBranchStatus("e_Py_pat",1);// e_Py_pat, &b_e_Py_pat);
   oldtree->SetBranchStatus("e_Pz_pat",1);// e_Pz_pat, &b_e_Pz_pat);
   oldtree->SetBranchStatus("e_En_pat",1);// e_En_pat, &b_e_En_pat);
   oldtree->SetBranchStatus("e_Eta_pat",1);// e_Eta_pat, &b_e_Eta_pat);
   oldtree->SetBranchStatus("e_Rapidity_pat",1);// e_Rapidity_pat, &b_e_Rapidity_pat);
   oldtree->SetBranchStatus("e_Phi_pat",1);// e_Phi_pat, &b_e_Phi_pat);
   oldtree->SetBranchStatus("e_Q_pat",1);// e_Q_pat, &b_e_Q_pat);
   oldtree->SetBranchStatus("e_EnSC_pat",1);// e_EnSC_pat, &b_e_EnSC_pat);
   oldtree->SetBranchStatus("e_EtSC_pat",1);// e_EtSC_pat, &b_e_EtSC_pat);
 //  oldtree->SetBranchStatus("e_dRbf_pat",1);// e_dRbf_pat, &b_e_dRbf_pat);
  // oldtree->SetBranchStatus("e_dR_pat",1);// e_dR_pat, &b_e_dR_pat);
   oldtree->SetBranchStatus("e_isoEcal_pat",1);// e_isoEcal_pat, &b_e_isoEcal_pat);
   oldtree->SetBranchStatus("e_isoHcal_pat",1);// e_isoHcal_pat, &b_e_isoHcal_pat);
   oldtree->SetBranchStatus("e_isoTrk_pat",1);// e_isoTrk_pat, &b_e_isoTrk_pat);
   oldtree->SetBranchStatus("e_isoComb_pat",1);// e_isoComb_pat, &b_e_isoComb_pat);
   oldtree->SetBranchStatus("e_Sigee_pat",1);// e_Sigee_pat, &b_e_Sigee_pat);
   oldtree->SetBranchStatus("e_SigIeIe_pat",1);// e_SigIeIe_pat, &b_e_SigIeIe_pat);
   oldtree->SetBranchStatus("e_dEtaIn_pat",1);// e_dEtaIn_pat, &b_e_dEtaIn_pat);
   oldtree->SetBranchStatus("e_dPhiIn_pat",1);// e_dPhiIn_pat, &b_e_dPhiIn_pat);
   oldtree->SetBranchStatus("e_HoE_pat",1);// e_HoE_pat, &b_e_HoE_pat);
   oldtree->SetBranchStatus("e_idLoose_pat",1);// e_idLoose_pat, &b_e_idLoose_pat);
   oldtree->SetBranchStatus("e_idTight_pat",1);// e_idTight_pat, &b_e_idTight_pat);
   oldtree->SetBranchStatus("e_idRobustHighEnergy_pat",1);// e_idRobustHighEnergy_pat, &b_e_idRobustHighEnergy_pat);
   oldtree->SetBranchStatus("e_idRobustLoose_pat",1);// e_idRobustLoose_pat, &b_e_idRobustLoose_pat);
   oldtree->SetBranchStatus("e_idRobustTight_pat",1);// e_idRobustTight_pat, &b_e_idRobustTight_pat);

   oldtree->SetBranchStatus("nJet_pat",1);// &nJet_pat, &b_nJet_pat);
   //oldtree->SetBranchStatus("jet_TMbf_pat",1);// jet_TMbf_pat, &b_jet_TMbf_pat);
  // oldtree->SetBranchStatus("jet_TM_pat",1);// jet_TM_pat, &b_jet_TM_pat);
  // oldtree->SetBranchStatus("jet_GJbf_pat",1);// jet_GJbf_pat, &b_jet_GJbf_pat);
   //oldtree->SetBranchStatus("jet_GJ_pat",1);// jet_GJ_pat, &b_jet_GJ_pat);
  // oldtree->SetBranchStatus("jet_RM_pat",1);// jet_RM_pat, &b_jet_RM_pat);
   oldtree->SetBranchStatus("jet_P_pat",1);// jet_P_pat, &b_jet_P_pat);
   oldtree->SetBranchStatus("jet_Pt_pat",1);// jet_Pt_pat, &b_jet_Pt_pat);
   oldtree->SetBranchStatus("jet_Px_pat",1);// jet_Px_pat, &b_jet_Px_pat);
   oldtree->SetBranchStatus("jet_Py_pat",1);// jet_Py_pat, &b_jet_Py_pat);
   oldtree->SetBranchStatus("jet_Pz_pat",1);// jet_Pz_pat, &b_jet_Pz_pat);
   oldtree->SetBranchStatus("jet_En_pat",1);// jet_En_pat, &b_jet_En_pat);
   oldtree->SetBranchStatus("jet_corrP_pat",1);// jet_corrP_pat, &b_jet_corrP_pat);
   oldtree->SetBranchStatus("jet_corrPt_pat",1);// jet_corrPt_pat, &b_jet_corrPt_pat);
   oldtree->SetBranchStatus("jet_corrPx_pat",1);// jet_corrPx_pat, &b_jet_corrPx_pat);
   oldtree->SetBranchStatus("jet_corrPy_pat",1);// jet_corrPy_pat, &b_jet_corrPy_pat);
   oldtree->SetBranchStatus("jet_corrPz_pat",1);// jet_corrPz_pat, &b_jet_corrPz_pat);
   oldtree->SetBranchStatus("jet_corrEn_pat",1);// jet_corrEn_pat, &b_jet_corrEn_pat);
   oldtree->SetBranchStatus("jet_Eta_pat",1);// jet_Eta_pat, &b_jet_Eta_pat);
   oldtree->SetBranchStatus("jet_Rapidity_pat",1);// jet_Rapidity_pat, &b_jet_Rapidity_pat);
   oldtree->SetBranchStatus("jet_Phi_pat",1);// jet_Phi_pat, &b_jet_Phi_pat);
   oldtree->SetBranchStatus("jet_Q_pat",1);// jet_Q_pat, &b_jet_Q_pat);
   oldtree->SetBranchStatus("jet_nTrk_pat",1);// jet_nTrk_pat, &b_jet_nTrk_pat);
   //oldtree->SetBranchStatus("jet_CaloPF_pat",1);// jet_CaloPF_pat, &b_jet_CaloPF_pat);
  // oldtree->SetBranchStatus("jet_dRbf_pat",1);// jet_dRbf_pat, &b_jet_dRbf_pat);
  // oldtree->SetBranchStatus("jet_dR_pat",1);// jet_dR_pat, &b_jet_dR_pat);
   //oldtree->SetBranchStatus("jet_GJdRbf_pat",1);// jet_GJdRbf_pat, &b_jet_GJdRbf_pat);
  // oldtree->SetBranchStatus("jet_GJdR_pat",1);// jet_GJdR_pat, &b_jet_GJdR_pat);
  // oldtree->SetBranchStatus("jet_RdR_pat",1);// jet_RdR_pat, &b_jet_RdR_pat);
   oldtree->SetBranchStatus("jet_cHadEn_pat",1);// jet_cHadEn_pat, &b_jet_cHadEn_pat);
   oldtree->SetBranchStatus("jet_cHadEF_pat",1);// jet_cHadEF_pat, &b_jet_cHadEF_pat);
   oldtree->SetBranchStatus("jet_nHadEn_pat",1);// jet_nHadEn_pat, &b_jet_nHadEn_pat);
   oldtree->SetBranchStatus("jet_nHadEF_pat",1);// jet_nHadEF_pat, &b_jet_nHadEF_pat);
   oldtree->SetBranchStatus("jet_cEmEn_pat",1);// jet_cEmEn_pat, &b_jet_cEmEn_pat);
   oldtree->SetBranchStatus("jet_cEmEF_pat",1);// jet_cEmEF_pat, &b_jet_cEmEF_pat);
   oldtree->SetBranchStatus("jet_nEmEn_pat",1);// jet_nEmEn_pat, &b_jet_nEmEn_pat);
   oldtree->SetBranchStatus("jet_nEmEF_pat",1);// jet_nEmEF_pat, &b_jet_nEmEF_pat);
   oldtree->SetBranchStatus("jet_muEn_pat",1);// jet_muEn_pat, &b_jet_muEn_pat);
   oldtree->SetBranchStatus("jet_muEF_pat",1);// jet_muEF_pat, &b_jet_muEF_pat);
   oldtree->SetBranchStatus("jet_cMult_pat",1);// jet_cMult_pat, &b_jet_cMult_pat);
   oldtree->SetBranchStatus("jet_nMult_pat",1);// jet_nMult_pat, &b_jet_nMult_pat);
   oldtree->SetBranchStatus("jet_muMult_pat",1);// jet_muMult_pat, &b_jet_muMult_pat);
   oldtree->SetBranchStatus("jet_n60_pat",1);// jet_n60_pat, &b_jet_n60_pat);
   oldtree->SetBranchStatus("jet_n90_pat",1);// jet_n90_pat, &b_jet_n90_pat);
   oldtree->SetBranchStatus("jet_totEM_pat",1);// jet_totEM_pat, &b_jet_totEM_pat);
   oldtree->SetBranchStatus("jet_totHad_pat",1);// jet_totHad_pat, &b_jet_totHad_pat);
   oldtree->SetBranchStatus("jet_EM_frac_pat",1);// jet_EM_frac_pat, &b_jet_EM_frac_pat);
   oldtree->SetBranchStatus("jet_EM_EB_frac_pat",1);// jet_EM_EB_frac_pat, &b_jet_EM_EB_frac_pat);
   oldtree->SetBranchStatus("jet_EM_EE_frac_pat",1);// jet_EM_EE_frac_pat, &b_jet_EM_EE_frac_pat);
   oldtree->SetBranchStatus("jet_EM_HF_frac_pat",1);// jet_EM_HF_frac_pat, &b_jet_EM_HF_frac_pat);
   oldtree->SetBranchStatus("jet_Had_HB_frac_pat",1);// jet_Had_HB_frac_pat, &b_jet_Had_HB_frac_pat);
   oldtree->SetBranchStatus("jet_Had_HE_frac_pat",1);// jet_Had_HE_frac_pat, &b_jet_Had_HE_frac_pat);
   oldtree->SetBranchStatus("jet_Had_HF_frac_pat",1);// jet_Had_HF_frac_pat, &b_jet_Had_HF_frac_pat);
   oldtree->SetBranchStatus("jet_Had_HO_frac_pat",1);// jet_Had_HO_frac_pat, &b_jet_Had_HO_frac_pat);



   oldtree->SetBranchStatus("nZ_pat",1);// &nZ_pat, &b_nZ_pat);
   oldtree->SetBranchStatus("z_dg1_pat",1);// z_dg1_pat, &b_z_dg1_pat);
   oldtree->SetBranchStatus("z_dg2_pat",1);// z_dg2_pat, &b_z_dg2_pat);
   oldtree->SetBranchStatus("z_M_pat",1);// z_M_pat, &b_z_M_pat);
   oldtree->SetBranchStatus("z_P_pat",1);// z_P_pat, &b_z_P_pat);
   oldtree->SetBranchStatus("z_Pt_pat",1);// z_Pt_pat, &b_z_Pt_pat);
   oldtree->SetBranchStatus("z_Px_pat",1);// z_Px_pat, &b_z_Px_pat);
   oldtree->SetBranchStatus("z_Py_pat",1);// z_Py_pat, &b_z_Py_pat);
   oldtree->SetBranchStatus("z_Pz_pat",1);// z_Pz_pat, &b_z_Pz_pat);
   oldtree->SetBranchStatus("z_En_pat",1);// z_En_pat, &b_z_En_pat);
   oldtree->SetBranchStatus("z_Eta_pat",1);// z_Eta_pat, &b_z_Eta_pat);
   oldtree->SetBranchStatus("z_Rapidity_pat",1);// z_Rapidity_pat, &b_z_Rapidity_pat);
   oldtree->SetBranchStatus("z_Phi_pat",1);// z_Phi_pat, &b_z_Phi_pat);
   oldtree->SetBranchStatus("z_dk_pat",1);// z_dk_pat, &b_z_dk_pat);

   TFile *newfile = new TFile("ttbar.root", "recreate");

   TTree *newTree = oldtree->CopyTree("nE_pat>1");

   TBranch *X = newTree->Branch("Xsection", &Xsection,"Xsection/D");
   TBranch *FE =newTree->Branch("filterEff", &FilterEff,"filterEff/D");
   TBranch *W = newTree->Branch("weight", &scale,"weight/D");
   cout<<scale<<endl;

   Int_t nEntries = (Int_t)newTree->GetEntries();
   cout<<"left entries:  "<<nEntries<<" ,   scale:      "<<scale<<endl;
   for (Int_t j=0; j<nEntries; j++) {
       X->Fill();
       FE->Fill();
       W->Fill();
   }
   
   long long max_tree_size = 200000000000LL; // 200 GB
   cout<<newTree->GetMaxTreeSize()<<endl;
   if(newTree->GetMaxTreeSize() < max_tree_size ) {
     newTree->SetMaxTreeSize(max_tree_size);
   }

   newfile->cd();
   newTree->Write();
   newfile->Close();
}
