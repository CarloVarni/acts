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

void ActsExamples::PersonalPlotTool::book(PersonalPlotTool::PersonalPlotCache& effPlotCache) const 
{
  PlotHelpers::Binning bPhi = m_cfg.varBinning.at("Phi");
  PlotHelpers::Binning bEta = m_cfg.varBinning.at("Eta");
  PlotHelpers::Binning bPt = m_cfg.varBinning.at("Pt");
  PlotHelpers::Binning bNRecoTimes = m_cfg.varBinning.at("N_Reco_Times");
  PlotHelpers::Binning bNContributingParticles = m_cfg.varBinning.at("N_Contributing_Particles");

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
  // Fakes
  effPlotCache.muon_fake_pT_vs_eta = PlotHelpers::bookHisto("muon_fake_pT_vs_eta",
							    "Fake Particle Pt vs #eta;Fake Particle Pt [Gev/C];Particle #eta",
							    bPt, bEta);
  effPlotCache.muon_fake_pT_vs_phi = PlotHelpers::bookHisto("muon_fake_pT_vs_phi",
							    "Fake Particle Pt vs #phi;Fake Particle Pt [GeV/C];Particle #phi",
							    bPt, bPhi);
  effPlotCache.muon_fake_eta_vs_phi = PlotHelpers::bookHisto("muon_fake_eta_vs_phi",
							     "Fake Particle #eta vs #phi;Particle #eta;Particle #phi",
							     bEta, bPhi);

  effPlotCache.muon_fake_pt = PlotHelpers::bookHisto("muon_fake_pt",
						     "Fake Particle Pt;Particle Pt [GeV/C];Entries",
						     bPt);
  effPlotCache.muon_fake_eta = PlotHelpers::bookHisto("muon_fake_eta",
						      "Fake Particle #eta;Particle #eta;Entries",
						      bEta);
  effPlotCache.muon_fake_phi = PlotHelpers::bookHisto("muon_fake_phi",
						      "Fake Particle #phi;Particle #phi;Entries",
						      bPhi);

  // Un-matched
  effPlotCache.muon_unmatched_pT_vs_eta = PlotHelpers::bookHisto("muon_unmatched_pT_vs_eta",
								 "Un-Matched Particle Pt vs #eta;Particle Pt [Gev/C];Particle #eta",
								 bPt, bEta);
  effPlotCache.muon_unmatched_pT_vs_phi = PlotHelpers::bookHisto("muon_unmatched_pT_vs_phi",
								 "Un-Matched Particle Pt vs #phi;Particle Pt [GeV/C];Particle #phi",
								 bPt, bPhi);
  effPlotCache.muon_unmatched_eta_vs_phi = PlotHelpers::bookHisto("muon_unmatched_eta_vs_phi",
								  "Un-Matched Particle #eta vs #phi;Particle #eta;Particle #phi",
								  bEta, bPhi);

  effPlotCache.muon_unmatched_pt = PlotHelpers::bookHisto("muon_unmatched_pt",
							  "Un-Matched Particle Pt;Particle Pt [GeV/C];Entries",
							  bPt);
  effPlotCache.muon_unmatched_eta = PlotHelpers::bookHisto("muon_unmatched_eta",
							   "Un-Matched Particle #eta;Particle #eta;Entries",
							   bEta);
  effPlotCache.muon_unmatched_phi = PlotHelpers::bookHisto("muon_unmatched_phi",
							   "Un-Matched Particle #phi;Particle #phi;Entries",
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
  // Other
  effPlotCache.n_reco_times = PlotHelpers::bookHisto("n_reco_times",
						     "Number of times was reconstructed;Entries",
						     bNRecoTimes);

  effPlotCache.n_contributing_particles = PlotHelpers::bookHisto("n_contributing_particles",
								 "Number of Contributing Particles;Entries",
								 bNContributingParticles);
}

void ActsExamples::PersonalPlotTool::clear(PersonalPlotCache& effPlotCache) const {

  delete effPlotCache.muon_fake_pT_vs_eta;
  delete effPlotCache.muon_fake_pT_vs_phi;
  delete effPlotCache.muon_fake_eta_vs_phi;

  delete effPlotCache.muon_fake_pt;
  delete effPlotCache.muon_fake_eta;
  delete effPlotCache.muon_fake_phi;


  delete effPlotCache.muon_unmatched_pT_vs_eta;
  delete effPlotCache.muon_unmatched_pT_vs_phi;
  delete effPlotCache.muon_unmatched_eta_vs_phi;

  delete effPlotCache.muon_unmatched_pt;
  delete effPlotCache.muon_unmatched_eta;
  delete effPlotCache.muon_unmatched_phi;



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
  delete effPlotCache.n_contributing_particles;
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


  effPlotCache.muon_fake_pT_vs_eta->Write();
  effPlotCache.muon_fake_pT_vs_phi->Write();
  effPlotCache.muon_fake_eta_vs_phi->Write();

  effPlotCache.muon_fake_pt->Write();
  effPlotCache.muon_fake_eta->Write();
  effPlotCache.muon_fake_phi->Write();

  effPlotCache.muon_unmatched_pT_vs_eta->Write();
  effPlotCache.muon_unmatched_pT_vs_phi->Write();
  effPlotCache.muon_unmatched_eta_vs_phi->Write();

  effPlotCache.muon_unmatched_pt->Write();
  effPlotCache.muon_unmatched_eta->Write();
  effPlotCache.muon_unmatched_phi->Write();


  effPlotCache.n_reco_times->Write();
  effPlotCache.n_contributing_particles->Write();
}

void ActsExamples::PersonalPlotTool::fill_nContributingParticles(PersonalPlotCache& effPlotCache,
								 std::size_t n_particles) const {
  PlotHelpers::fillHisto(effPlotCache.n_contributing_particles, n_particles); 
}

void ActsExamples::PersonalPlotTool::fill_fake(PersonalPlotCache& effPlotCache,
					       float t_pt, float t_eta, float t_phi) const {

  PlotHelpers::fillHisto(effPlotCache.muon_fake_pT_vs_eta , t_pt, t_eta);
  PlotHelpers::fillHisto(effPlotCache.muon_fake_pT_vs_phi , t_pt, t_phi);
  PlotHelpers::fillHisto(effPlotCache.muon_fake_eta_vs_phi, t_eta, t_phi);

  PlotHelpers::fillHisto(effPlotCache.muon_fake_pt, t_pt);
  PlotHelpers::fillHisto(effPlotCache.muon_fake_eta, t_eta);
  PlotHelpers::fillHisto(effPlotCache.muon_fake_phi, t_phi);
}

void ActsExamples::PersonalPlotTool::fill_unmatched(PersonalPlotCache& effPlotCache,
						    float t_pt, float t_eta, float t_phi) const {

  PlotHelpers::fillHisto(effPlotCache.muon_unmatched_pT_vs_eta, t_pt, t_eta);
  PlotHelpers::fillHisto(effPlotCache.muon_unmatched_pT_vs_phi, t_pt, t_phi);
  PlotHelpers::fillHisto(effPlotCache.muon_unmatched_eta_vs_phi, t_eta, t_phi);

  PlotHelpers::fillHisto(effPlotCache.muon_unmatched_pt, t_pt);
  PlotHelpers::fillHisto(effPlotCache.muon_unmatched_eta, t_eta);
  PlotHelpers::fillHisto(effPlotCache.muon_unmatched_phi, t_phi);
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
