// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Validation/PersonalPlotTool.hpp"

#include "Acts/Utilities/Helpers.hpp"

using Acts::VectorHelpers::eta;
using Acts::VectorHelpers::perp;
using Acts::VectorHelpers::phi;

ActsExamples::PersonalPlotTool::PersonalPlotTool(
						 const ActsExamples::PersonalPlotTool::Config& cfg, Acts::Logging::Level lvl)
    : m_cfg(cfg), m_logger(Acts::getDefaultLogger("PersonalPlotTool", lvl)) {}

void ActsExamples::PersonalPlotTool::book(
    PersonalPlotTool::PersonalPlotCache& effPlotCache) const {
  PlotHelpers::Binning bPhi = m_cfg.varBinning.at("Phi");
  PlotHelpers::Binning bEta = m_cfg.varBinning.at("Eta");
  PlotHelpers::Binning bPt = m_cfg.varBinning.at("Pt");
  PlotHelpers::Binning bNRecoTimes = m_cfg.varBinning.at("N_Reco_Times");

  ACTS_DEBUG("Initialize the histograms for efficiency plots");
  effPlotCache.trackEff_vs_pT_vs_eta =
    PlotHelpers::bookEff("eff_muon_pt_vs_eta",
			 "Tracking Efficiency;Truth pT [GeV/c];Truth #eta",
			 bPt, bEta);
  
  effPlotCache.trackEff_vs_pT_vs_phi =
    PlotHelpers::bookEff("eff_muon_pt_vs_phi",
			 "Tracking Efficiency;Truth pT [GeV/c];Truth #phi",
			 bPt, bPhi);
  
  effPlotCache.trackEff_vs_eta_vs_phi =
    PlotHelpers::bookEff("eff_muon_eta_vs_phi",
			 "Tracking Efficiency;Truth #eta;Truth #phi",
			 bEta, bPhi);
  
  // Reco
  effPlotCache.muon_reco_pT_vs_eta =
    PlotHelpers::bookHisto("muon_reco_pT_vs_eta",
			   "Reco Muon pT [GeV/c];Reco Muon #eta",
			   bPt, bEta);

  effPlotCache.muon_reco_pT_vs_phi =
    PlotHelpers::bookHisto("muon_reco_pT_vs_phi",
			   "Reco Muon pT [GeV/c];Reco Muon #phi",
			   bPt, bPhi);

  effPlotCache.muon_reco_eta_vs_phi =
    PlotHelpers::bookHisto("muon_reco_eta_vs_phi",
			   "Reco Muon #eta;Reco Muon #phi",
			   bEta, bPhi);

  effPlotCache.muon_reco_pt = PlotHelpers::bookHisto("muon_reco_pt", 
						     "Reco Muon pT [GeV/c];Entries",
						     bPt);
  effPlotCache.muon_reco_eta = PlotHelpers::bookHisto("muon_reco_eta", 
						      "Reco Muon #eta [GeV/c];Entries",
						      bEta);
  effPlotCache.muon_reco_phi = PlotHelpers::bookHisto("muon_reco_phi", 
						      "Reco Muon #phi [GeV/c];Entries",
						      bPhi);

  // Truth  
  effPlotCache.muon_truth_pT_vs_eta =
    PlotHelpers::bookHisto("muon_truth_pT_vs_eta",
			   "Truth Muon pT [GeV/c];Truth Muon #eta",
			   bPt, bEta);

  effPlotCache.muon_truth_pT_vs_phi =
    PlotHelpers::bookHisto("muon_truth_pT_vs_phi",
			   "Truth Muon pT [GeV/c];Truth Muon #phi",
			   bPt, bPhi);

  effPlotCache.muon_truth_eta_vs_phi =
    PlotHelpers::bookHisto("muon_truth_eta_vs_phi",
			   "Truth Muon #eta;Truth Muon #phi",
			   bEta, bPhi);

  effPlotCache.muon_truth_pt = PlotHelpers::bookHisto("muon_truth_pt", 
						      "Truth Muon pT [GeV/c];Entries",
						      bPt);
  effPlotCache.muon_truth_eta = PlotHelpers::bookHisto("muon_truth_eta", 
						       "Truth Muon #eta [GeV/c];Entries",
						       bEta);
  effPlotCache.muon_truth_phi = PlotHelpers::bookHisto("muon_truth_phi", 
						       "Truth Muon #phi [GeV/c];Entries",
						       bPhi);

  effPlotCache.n_reco_times = PlotHelpers::bookHisto("n_reco_times",
						     "Number of times was reconstructed;Entries",
						     bNRecoTimes);

}

void ActsExamples::PersonalPlotTool::clear(PersonalPlotCache& effPlotCache) const {
  delete effPlotCache.trackEff_vs_pT_vs_eta;
  delete effPlotCache.trackEff_vs_pT_vs_phi;
  delete effPlotCache.trackEff_vs_eta_vs_phi;

  delete effPlotCache.muon_reco_pT_vs_eta;
  delete effPlotCache.muon_reco_pT_vs_phi;
  delete effPlotCache.muon_reco_eta_vs_phi;

  delete effPlotCache.muon_reco_pt;
  delete effPlotCache.muon_reco_eta;
  delete effPlotCache.muon_reco_phi;

  delete effPlotCache.muon_truth_pT_vs_eta;
  delete effPlotCache.muon_truth_pT_vs_phi;
  delete effPlotCache.muon_truth_eta_vs_phi;

  delete effPlotCache.muon_truth_pt;
  delete effPlotCache.muon_truth_eta;
  delete effPlotCache.muon_truth_phi;

  delete effPlotCache.n_reco_times;
}

void ActsExamples::PersonalPlotTool::write(
    const PersonalPlotTool::PersonalPlotCache& effPlotCache) const {
  ACTS_DEBUG("Write the plots to output file.");
  effPlotCache.trackEff_vs_pT_vs_eta->Write();
  effPlotCache.trackEff_vs_pT_vs_phi->Write();
  effPlotCache.trackEff_vs_eta_vs_phi->Write();

  effPlotCache.muon_reco_pT_vs_eta->Write();
  effPlotCache.muon_reco_pT_vs_phi->Write();
  effPlotCache.muon_reco_eta_vs_phi->Write();

  effPlotCache.muon_reco_pt->Write();
  effPlotCache.muon_reco_eta->Write();
  effPlotCache.muon_reco_phi->Write();

  effPlotCache.muon_truth_pT_vs_eta->Write();
  effPlotCache.muon_truth_pT_vs_phi->Write();
  effPlotCache.muon_truth_eta_vs_phi->Write();

  effPlotCache.muon_truth_pt->Write();
  effPlotCache.muon_truth_eta->Write();
  effPlotCache.muon_truth_phi->Write();

  effPlotCache.n_reco_times->Write();
}

void ActsExamples::PersonalPlotTool::fill_reco(PersonalPlotCache& effPlotCache,
					       float t_pT, float t_eta, float t_phi) const {
  PlotHelpers::fillHisto(effPlotCache.muon_reco_pT_vs_eta, t_pT, t_eta);
  PlotHelpers::fillHisto(effPlotCache.muon_reco_pT_vs_phi, t_pT, t_phi);
  PlotHelpers::fillHisto(effPlotCache.muon_reco_eta_vs_phi, t_eta, t_phi);

  PlotHelpers::fillHisto(effPlotCache.muon_reco_pt, t_pT);
  PlotHelpers::fillHisto(effPlotCache.muon_reco_eta, t_eta);
  PlotHelpers::fillHisto(effPlotCache.muon_reco_phi, t_phi);

}

void ActsExamples::PersonalPlotTool::fill_truth(PersonalPlotTool::PersonalPlotCache& effPlotCache,
						const ActsFatras::Particle& truthParticle,
						int RecoTimes,
						bool status) const {

  const auto t_phi = phi(truthParticle.unitDirection());
  const auto t_eta = eta(truthParticle.unitDirection());
  const auto t_pT = truthParticle.transverseMomentum();
  
  PlotHelpers::fillEff(effPlotCache.trackEff_vs_pT_vs_eta, t_pT, t_eta, status );
  PlotHelpers::fillEff(effPlotCache.trackEff_vs_pT_vs_phi, t_pT, t_phi, status );
  PlotHelpers::fillEff(effPlotCache.trackEff_vs_eta_vs_phi, t_eta, t_phi, status );
  
  PlotHelpers::fillHisto(effPlotCache.muon_truth_pT_vs_eta, t_pT, t_eta);
  PlotHelpers::fillHisto(effPlotCache.muon_truth_pT_vs_phi, t_pT, t_phi);
  PlotHelpers::fillHisto(effPlotCache.muon_truth_eta_vs_phi, t_eta, t_phi);

  PlotHelpers::fillHisto(effPlotCache.muon_truth_pt, t_pT);
  PlotHelpers::fillHisto(effPlotCache.muon_truth_eta, t_eta);
  PlotHelpers::fillHisto(effPlotCache.muon_truth_phi, t_phi);

  PlotHelpers::fillHisto(effPlotCache.n_reco_times, RecoTimes);
}
