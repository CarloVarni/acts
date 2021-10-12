// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Validation/PhysicsPerformancePlotTool.hpp"

#include "Acts/Utilities/Helpers.hpp"

using Acts::VectorHelpers::eta;
using Acts::VectorHelpers::perp;
using Acts::VectorHelpers::phi;

namespace ActsExamples {

  PhysicsPerformancePlotTool::PhysicsPerformancePlotTool(const PhysicsPerformancePlotTool::Config& cfg,
							 Acts::Logging::Level lvl)
    : m_cfg(cfg), 
      m_logger(Acts::getDefaultLogger(m_cfg.tool_name.c_str(), lvl)) 
  {
    if (m_cfg.tool_name.empty())
      throw std::invalid_argument("Missing tool name");
    if (m_cfg.particle_name.empty())
      throw std::invalid_argument("Missing particle name");
  }

  void 
  PhysicsPerformancePlotTool::book(PhysicsPerformancePlotTool::PhysicsPerformancePlotCache& effPlotCache) const 
  {
    std::string particle_name = m_cfg.particle_name;
    PlotHelpers::Binning bPhi = m_cfg.varBinning.at("Phi");
    PlotHelpers::Binning bEta = m_cfg.varBinning.at("Eta");
    PlotHelpers::Binning bPt = m_cfg.varBinning.at("Pt");
    PlotHelpers::Binning bErrPt = m_cfg.varBinning.at("err_Pt");
    PlotHelpers::Binning bD0 = m_cfg.varBinning.at("d0");
    PlotHelpers::Binning bZ0 = m_cfg.varBinning.at("z0");
    PlotHelpers::Binning bErrD0 = m_cfg.varBinning.at("err_d0");
    PlotHelpers::Binning bErrZ0 = m_cfg.varBinning.at("err_z0");
    PlotHelpers::Binning bSignificanceD0 = m_cfg.varBinning.at("significance_d0");
    PlotHelpers::Binning bSignificanceZ0 = m_cfg.varBinning.at("significance_z0");
    PlotHelpers::Binning bNRecoTimes = m_cfg.varBinning.at("N_Reco_Times");
    PlotHelpers::Binning bNContributingParticles = m_cfg.varBinning.at("N_Contributing_Particles");
    
    ACTS_DEBUG("Initialize the histograms for efficiency plots");
    effPlotCache.trackEff_vs_pT_vs_eta =
      PlotHelpers::bookEff((particle_name + "_eff_pt_vs_eta").c_str(),
			   (particle_name + "Tracking Efficiency;Truth pT [GeV/c];Truth #eta").c_str(),
			   bPt, bEta);
    
    effPlotCache.trackEff_vs_pT_vs_phi =
      PlotHelpers::bookEff((particle_name + "_eff_pt_vs_phi").c_str(),
			   (particle_name + "Tracking Efficiency;Truth pT [GeV/c];Truth #phi").c_str(),
			   bPt, bPhi);
    
    effPlotCache.trackEff_vs_eta_vs_phi =
      PlotHelpers::bookEff((particle_name + "_eff_eta_vs_phi").c_str(),
			   (particle_name + "Tracking Efficiency;Truth #eta;Truth #phi").c_str(),
			   bEta, bPhi);
    
    effPlotCache.particle_pT_vs_eta =
      PlotHelpers::bookHisto((particle_name + "_pT_vs_eta").c_str(),
			     (particle_name + "Particle pT [GeV/c];Particle #eta").c_str(),
			     bPt, bEta);
    
    effPlotCache.particle_pT_vs_phi =
      PlotHelpers::bookHisto((particle_name + "_pT_vs_phi").c_str(),
			     "Particle pT [GeV/c];Particle #phi",
			     bPt, bPhi);
    
    effPlotCache.particle_eta_vs_phi =
      PlotHelpers::bookHisto((particle_name + "_eta_vs_phi").c_str(),
			     "Particle #eta;Particle #phi",
			     bEta, bPhi);
    
    effPlotCache.particle_pt = PlotHelpers::bookHisto((particle_name + "_pt").c_str(), 
						      "Particle pT [GeV/c];Entries",
						      bPt);
    effPlotCache.particle_eta = PlotHelpers::bookHisto((particle_name + "_eta").c_str(), 
						       "Particle #eta [GeV/c];Entries",
						       bEta);
    effPlotCache.particle_phi = PlotHelpers::bookHisto((particle_name + "_phi").c_str(), 
						       "Particle #phi [GeV/c];Entries",
						       bPhi);
 
    effPlotCache.profile_particle_pT_vs_eta = PlotHelpers::bookProf((particle_name + "_profile_pt_vs_eta").c_str(),
								    "Particle Pt [GeV/C]",
								    bEta, bPt);
    effPlotCache.profile_particle_pT_vs_phi = PlotHelpers::bookProf((particle_name + "_profile_pt_vs_phi").c_str(),
								    "Particle Pt [GeV/C]",
								    bPhi, bPt);

    effPlotCache.profile_particle_pTqOverP_vs_eta = PlotHelpers::bookProf((particle_name + "_profile_pTqOverP_vs_eta").c_str(),
									  "Particle Pt x sigma(Q over P)",
									  bEta, bPt);
    effPlotCache.profile_particle_pTqOverP_vs_phi = PlotHelpers::bookProf((particle_name + "_profile_pTqOverP_vs_phi").c_str(),
									  "Particle Pt x sigma(Q over P)",
									  bPhi, bPt);

    effPlotCache.profile_particle_err_inv_pT_vs_eta = PlotHelpers::bookProf((particle_name + "_profile_err_inv_pt_vs_eta").c_str(),
									    "Particle sigma inv Pt ",
									    bEta, bErrPt);

    effPlotCache.profile_particle_err_inv_pT_vs_phi = PlotHelpers::bookProf((particle_name + "_profile_err_inv_pt_vs_phi").c_str(),
									    "Particle sigma inv Pt ",
									    bPhi, bErrPt);

    //
    effPlotCache.particle_d0 = PlotHelpers::bookHisto((particle_name + "_d0").c_str(),
						      "d0 [mm]; Entries",
						      bD0);
    effPlotCache.particle_z0 = PlotHelpers::bookHisto((particle_name + "_z0").c_str(),
						      "z0 [mm]; Entries",
						      bZ0);
    effPlotCache.particle_err_d0 = PlotHelpers::bookHisto((particle_name + "_err_d0").c_str(),
							  "#sigma d0 [mm]; Entries",
							  bErrD0);
    effPlotCache.particle_err_z0 = PlotHelpers::bookHisto((particle_name + "_err_z0").c_str(),
							  "#sigma z0 [mm]; Entries",
							  bErrZ0);
    effPlotCache.particle_significance_d0 = PlotHelpers::bookHisto((particle_name + "_significance_d0").c_str(),
								   "Significance d0; Entries",
								   bSignificanceD0);
    effPlotCache.particle_significance_z0 = PlotHelpers::bookHisto((particle_name + "_significance_z0").c_str(),
								   "Significance z0; Entries",
								   bSignificanceZ0);
    

    effPlotCache.profile_particle_d0_vs_eta = PlotHelpers::bookProf((particle_name + "_profile_d0_vs_eta").c_str(),
								    "Impact Parameter d0",
								    bEta, bD0);
    effPlotCache.profile_particle_d0_vs_phi = PlotHelpers::bookProf((particle_name + "_profile_d0_vs_phi").c_str(),
								    "Impact Parameter d0",
								    bPhi, bD0);
    effPlotCache.profile_particle_z0_vs_eta = PlotHelpers::bookProf((particle_name + "_profile_z0_vs_eta").c_str(),
								    "Impact Parameter z0",
								    bEta, bZ0);
    effPlotCache.profile_particle_z0_vs_phi = PlotHelpers::bookProf((particle_name + "_profile_z0_vs_phi").c_str(),
								    "Impact Parameter z0",
								    bPhi, bZ0);

    effPlotCache.profile_particle_err_d0_vs_eta = PlotHelpers::bookProf((particle_name + "_profile_err_d0_vs_eta").c_str(),
									"Impact Parameter Error d0",
									bEta, bErrD0);
    effPlotCache.profile_particle_err_d0_vs_phi = PlotHelpers::bookProf((particle_name + "_profile_err_d0_vs_phi").c_str(),
									"Impact Parameter Error d0",
									bPhi, bErrD0);
    effPlotCache.profile_particle_err_z0_vs_eta = PlotHelpers::bookProf((particle_name + "_profile_err_z0_vs_eta").c_str(),
									"Impact Parameter Error z0",
									bEta, bErrZ0);
    effPlotCache.profile_particle_err_z0_vs_phi = PlotHelpers::bookProf((particle_name + "_profile_err_z0_vs_phi").c_str(),
									"Impact Parameter Error z0",
									bPhi, bErrZ0);
    
    // 2D
    effPlotCache.particle_d0_vs_eta = PlotHelpers::bookHisto((particle_name + "_d0_vs_eta").c_str(),
							     "d0 [mm]; #eta",
							     bD0, bEta);
    effPlotCache.particle_d0_vs_phi = PlotHelpers::bookHisto((particle_name + "_d0_vs_phi").c_str(),
							     "d0 [mm]; #phi",
							     bD0, bPhi);
    effPlotCache.particle_z0_vs_eta = PlotHelpers::bookHisto((particle_name + "_z0_vs_eta").c_str(),
							     "z0 [mm]; #eta",
							     bZ0, bEta);
    effPlotCache.particle_z0_vs_phi = PlotHelpers::bookHisto((particle_name + "_z0_vs_phi").c_str(),
							     "z0 [mm]; #phi",
							     bZ0, bPhi);

    effPlotCache.particle_err_d0_vs_eta = PlotHelpers::bookHisto((particle_name + "_err_d0_vs_eta").c_str(),
								 "#sigma d0 [mm]; #eta",
								 bErrD0, bEta);
    effPlotCache.particle_err_d0_vs_phi = PlotHelpers::bookHisto((particle_name + "_err_d0_vs_phi").c_str(),
								 "#sigma d0 [mm]; #phi",
								 bErrD0, bPhi);
    effPlotCache.particle_err_z0_vs_eta = PlotHelpers::bookHisto((particle_name + "_err_z0_vs_eta").c_str(),
								 "#sigma z0 [mm]; #eta",
								 bErrZ0, bEta);
    effPlotCache.particle_err_z0_vs_phi = PlotHelpers::bookHisto((particle_name + "_err_z0_vs_phi").c_str(),
								 "#sigma z0 [mm]; #phi",
								 bErrZ0, bPhi);

    effPlotCache.particle_significance_d0_vs_eta = PlotHelpers::bookHisto((particle_name + "_significance_d0_vs_eta").c_str(),
									  "Significance d0; #eta",
									  bSignificanceD0, bEta);
    effPlotCache.particle_significance_d0_vs_phi = PlotHelpers::bookHisto((particle_name + "_significance_d0_vs_phi").c_str(),
									  "Significance d0; #phi",
									  bSignificanceD0, bPhi);
    effPlotCache.particle_significance_z0_vs_eta = PlotHelpers::bookHisto((particle_name + "_significance_z0_vs_eta").c_str(),
									  "Significance z0; #eta",
									  bSignificanceZ0, bEta);
    effPlotCache.particle_significance_z0_vs_phi = PlotHelpers::bookHisto((particle_name + "_significance_z0_vs_phi").c_str(),
									  "Significance z0; #phi",
									  bSignificanceZ0, bPhi);

    // Other
    effPlotCache.n_reco_times = PlotHelpers::bookHisto((particle_name + "_n_reco_times").c_str(),
						       "Number of times the particle was reconstructed;Entries",
						       bNRecoTimes);
    
    effPlotCache.n_contributing_particles = PlotHelpers::bookHisto((particle_name + "_n_contributing_particles").c_str(),
								   "Number of Contributing Particles;Entries",
								   bNContributingParticles);
  }
  
  void
  PhysicsPerformancePlotTool::clear(PhysicsPerformancePlotTool::PhysicsPerformancePlotCache& effPlotCache) const 
  {  
    delete effPlotCache.particle_pT_vs_eta;
    delete effPlotCache.particle_pT_vs_phi;
    delete effPlotCache.particle_eta_vs_phi;
    
    delete effPlotCache.particle_pt;
    delete effPlotCache.particle_eta;
    delete effPlotCache.particle_phi;
    delete effPlotCache.profile_particle_pT_vs_eta;
    delete effPlotCache.profile_particle_pT_vs_phi;
    
    delete effPlotCache.profile_particle_err_inv_pT_vs_eta;
    delete effPlotCache.profile_particle_err_inv_pT_vs_phi;

    delete effPlotCache.profile_particle_pTqOverP_vs_eta;
    delete effPlotCache.profile_particle_pTqOverP_vs_phi;

    delete effPlotCache.particle_d0;
    delete effPlotCache.particle_z0;
    delete effPlotCache.particle_err_d0;
    delete effPlotCache.particle_err_z0;
    delete effPlotCache.particle_significance_d0;
    delete effPlotCache.particle_significance_z0;

    delete effPlotCache.particle_d0_vs_eta;
    delete effPlotCache.particle_d0_vs_phi;
    delete effPlotCache.particle_z0_vs_eta;
    delete effPlotCache.particle_z0_vs_phi;

    delete effPlotCache.profile_particle_d0_vs_eta;
    delete effPlotCache.profile_particle_d0_vs_phi;
    delete effPlotCache.profile_particle_z0_vs_eta;
    delete effPlotCache.profile_particle_z0_vs_phi;

    delete effPlotCache.profile_particle_err_d0_vs_eta;
    delete effPlotCache.profile_particle_err_d0_vs_phi;
    delete effPlotCache.profile_particle_err_z0_vs_eta;
    delete effPlotCache.profile_particle_err_z0_vs_phi;

    delete effPlotCache.particle_err_d0_vs_eta;
    delete effPlotCache.particle_err_d0_vs_phi;
    delete effPlotCache.particle_err_z0_vs_eta;
    delete effPlotCache.particle_err_z0_vs_phi;

    delete effPlotCache.particle_significance_d0_vs_eta;
    delete effPlotCache.particle_significance_d0_vs_phi;
    delete effPlotCache.particle_significance_z0_vs_eta;
    delete effPlotCache.particle_significance_z0_vs_phi;

    delete effPlotCache.n_reco_times;
    delete effPlotCache.n_contributing_particles;
  }
  
  void
  PhysicsPerformancePlotTool::write(const PhysicsPerformancePlotTool::PhysicsPerformancePlotCache& effPlotCache) const 
  {  
    ACTS_DEBUG("Write the plots to output file.");
    
    effPlotCache.trackEff_vs_pT_vs_eta->Write();
    effPlotCache.trackEff_vs_pT_vs_phi->Write();
    effPlotCache.trackEff_vs_eta_vs_phi->Write();
    
    effPlotCache.particle_pT_vs_eta->Write();
    effPlotCache.particle_pT_vs_phi->Write();
    effPlotCache.particle_eta_vs_phi->Write();
    
    effPlotCache.particle_pt->Write();
    effPlotCache.particle_eta->Write();
    effPlotCache.particle_phi->Write();
    effPlotCache.profile_particle_pT_vs_eta->Write();
    effPlotCache.profile_particle_pT_vs_phi->Write();
    
    effPlotCache.profile_particle_err_inv_pT_vs_eta->Write();
    effPlotCache.profile_particle_err_inv_pT_vs_phi->Write();

    effPlotCache.profile_particle_pTqOverP_vs_eta->Write();
    effPlotCache.profile_particle_pTqOverP_vs_phi->Write();

    effPlotCache.particle_d0->Write();
    effPlotCache.particle_z0->Write();
    effPlotCache.particle_err_d0->Write();
    effPlotCache.particle_err_z0->Write();
    effPlotCache.particle_significance_d0->Write();
    effPlotCache.particle_significance_z0->Write();

    effPlotCache.particle_d0_vs_eta->Write();
    effPlotCache.particle_d0_vs_phi->Write();
    effPlotCache.particle_z0_vs_eta->Write();
    effPlotCache.particle_z0_vs_phi->Write();

    effPlotCache.profile_particle_d0_vs_eta->Write();
    effPlotCache.profile_particle_d0_vs_phi->Write();
    effPlotCache.profile_particle_z0_vs_eta->Write();
    effPlotCache.profile_particle_z0_vs_phi->Write();

    effPlotCache.profile_particle_err_d0_vs_eta->Write();
    effPlotCache.profile_particle_err_d0_vs_phi->Write();
    effPlotCache.profile_particle_err_z0_vs_eta->Write();
    effPlotCache.profile_particle_err_z0_vs_phi->Write();

    effPlotCache.particle_err_d0_vs_eta->Write();
    effPlotCache.particle_err_d0_vs_phi->Write();
    effPlotCache.particle_err_z0_vs_eta->Write();
    effPlotCache.particle_err_z0_vs_phi->Write();

    effPlotCache.particle_significance_d0_vs_eta->Write();
    effPlotCache.particle_significance_d0_vs_phi->Write();
    effPlotCache.particle_significance_z0_vs_eta->Write();
    effPlotCache.particle_significance_z0_vs_phi->Write();

    effPlotCache.n_reco_times->Write();
    effPlotCache.n_contributing_particles->Write();
  }
  
  void
  PhysicsPerformancePlotTool::fill_nContributingParticles(PhysicsPerformancePlotTool::PhysicsPerformancePlotCache& effPlotCache,
							  std::size_t n_particles) const 
  {
    PlotHelpers::fillHisto(effPlotCache.n_contributing_particles, n_particles); 
  }

  void
  PhysicsPerformancePlotTool::fill(PhysicsPerformancePlotTool::PhysicsPerformancePlotCache& effPlotCache,
				   float t_pt, float t_eta, float t_phi) const 
  {
    PlotHelpers::fillHisto(effPlotCache.particle_pT_vs_eta , t_pt, t_eta);
    PlotHelpers::fillHisto(effPlotCache.particle_pT_vs_phi , t_pt, t_phi);
    PlotHelpers::fillHisto(effPlotCache.particle_eta_vs_phi, t_eta, t_phi);
    
    PlotHelpers::fillHisto(effPlotCache.particle_pt, t_pt);
    PlotHelpers::fillHisto(effPlotCache.particle_eta, t_eta);
    PlotHelpers::fillHisto(effPlotCache.particle_phi, t_phi);

    PlotHelpers::fillProf(effPlotCache.profile_particle_pT_vs_eta, t_eta, t_pt);
    PlotHelpers::fillProf(effPlotCache.profile_particle_pT_vs_phi, t_phi, t_pt);
  }

  void 
  PhysicsPerformancePlotTool::fill(PhysicsPerformancePlotCache& effPlotCache,
				   const TrackParameters& trackParams,
				   std::size_t nContributingParticles) const
  {
    this->fill(effPlotCache, trackParams);
    this->fill_nContributingParticles(effPlotCache, nContributingParticles);
  }
  
  void 
  PhysicsPerformancePlotTool::fill(PhysicsPerformancePlotCache& effPlotCache,
				   const TrackParameters& trackParams) const 
  {
    if (not trackParams.covariance().has_value())
      throw std::runtime_error("No covariance in track parameter object");

    const auto& momentum = trackParams.momentum();
    const auto pT = Acts::VectorHelpers::perp(momentum);
    const auto eta = Acts::VectorHelpers::eta(momentum);
    const auto phi = Acts::VectorHelpers::phi(momentum);

    this->fill(effPlotCache, pT, eta, phi);

    const auto& params = trackParams.parameters();
    const auto& cov_matrix = trackParams.covariance().value();

    // Acts::eBoundLoc0, Acts::eBoundLoc1, Acts::eBoundPhi, Acts::eBoundTheta, Acts::eBoundQOverP, Acts::eBoundTime
    const auto d0 = params[Acts::eBoundLoc0];
    const auto z0 = params[Acts::eBoundLoc1];
    const auto theta = params[Acts::eBoundTheta];
    const auto err_d0 = std::sqrt( cov_matrix(Acts::eBoundLoc0, Acts::eBoundLoc0) );
    const auto err_z0 = std::sqrt( cov_matrix(Acts::eBoundLoc1, Acts::eBoundLoc1) );
    const auto err_qOverP = sqrt( cov_matrix(Acts::eBoundQOverP, Acts::eBoundQOverP) ); 
    const auto err_pT = pT * pT * err_qOverP;

    PlotHelpers::fillProf(effPlotCache.profile_particle_pTqOverP_vs_eta, eta, err_qOverP);
    PlotHelpers::fillProf(effPlotCache.profile_particle_pTqOverP_vs_phi, phi, err_qOverP);

    PlotHelpers::fillProf(effPlotCache.profile_particle_err_inv_pT_vs_eta, eta, err_pT);
    PlotHelpers::fillProf(effPlotCache.profile_particle_err_inv_pT_vs_phi, phi, err_pT);

    PlotHelpers::fillHisto(effPlotCache.particle_d0, d0);
    PlotHelpers::fillHisto(effPlotCache.particle_z0, z0 * sin(theta));
    PlotHelpers::fillHisto(effPlotCache.particle_err_d0, err_d0);
    PlotHelpers::fillHisto(effPlotCache.particle_err_z0, err_z0 * sin(theta));
    PlotHelpers::fillHisto(effPlotCache.particle_significance_d0, err_d0 / d0 * 1.);
    PlotHelpers::fillHisto(effPlotCache.particle_significance_z0, err_z0 / z0 * 1.);

    PlotHelpers::fillHisto(effPlotCache.particle_d0_vs_eta, d0, eta);
    PlotHelpers::fillHisto(effPlotCache.particle_d0_vs_phi, d0, phi);
    PlotHelpers::fillHisto(effPlotCache.particle_z0_vs_eta, z0 * sin(theta), eta);
    PlotHelpers::fillHisto(effPlotCache.particle_z0_vs_phi, z0 * sin(theta), phi);

    PlotHelpers::fillProf(effPlotCache.profile_particle_d0_vs_eta, eta, d0);
    PlotHelpers::fillProf(effPlotCache.profile_particle_d0_vs_phi, phi, d0);
    PlotHelpers::fillProf(effPlotCache.profile_particle_z0_vs_eta, eta, z0 * sin(theta));
    PlotHelpers::fillProf(effPlotCache.profile_particle_z0_vs_phi, phi, z0 * sin(theta));

    PlotHelpers::fillProf(effPlotCache.profile_particle_err_d0_vs_eta, eta, err_d0);
    PlotHelpers::fillProf(effPlotCache.profile_particle_err_d0_vs_phi, phi, err_d0);
    PlotHelpers::fillProf(effPlotCache.profile_particle_err_z0_vs_eta, eta, err_z0 * sin(theta));
    PlotHelpers::fillProf(effPlotCache.profile_particle_err_z0_vs_phi, phi, err_z0 * sin(theta));

    PlotHelpers::fillHisto(effPlotCache.particle_err_d0_vs_eta, err_d0, eta);
    PlotHelpers::fillHisto(effPlotCache.particle_err_d0_vs_phi, err_d0, phi);
    PlotHelpers::fillHisto(effPlotCache.particle_err_z0_vs_eta, err_z0 * sin(theta), eta);
    PlotHelpers::fillHisto(effPlotCache.particle_err_z0_vs_phi, err_z0 * sin(theta), phi);

    PlotHelpers::fillHisto(effPlotCache.particle_significance_d0_vs_eta, err_d0 / d0 * 1., eta);
    PlotHelpers::fillHisto(effPlotCache.particle_significance_d0_vs_phi, err_d0 / d0 * 1., phi);
    PlotHelpers::fillHisto(effPlotCache.particle_significance_z0_vs_eta, err_z0 / z0 * 1., eta);
    PlotHelpers::fillHisto(effPlotCache.particle_significance_z0_vs_phi, err_z0 / z0 * 1., phi);
  }

  void
  PhysicsPerformancePlotTool::fill(PhysicsPerformancePlotTool::PhysicsPerformancePlotCache& effPlotCache,
				   const ActsFatras::Particle& truthParticle,
				   int RecoTimes,
				   bool status) const 
  {
    const auto t_pT = truthParticle.transverseMomentum();
    const auto t_eta = eta(truthParticle.unitDirection());
    const auto t_phi = phi(truthParticle.unitDirection());

    PlotHelpers::fillEff(effPlotCache.trackEff_vs_pT_vs_eta, t_pT, t_eta, status );
    PlotHelpers::fillEff(effPlotCache.trackEff_vs_pT_vs_phi, t_pT, t_phi, status );
    PlotHelpers::fillEff(effPlotCache.trackEff_vs_eta_vs_phi, t_eta, t_phi, status );

    this->fill(effPlotCache, t_pT, t_eta, t_phi);
    PlotHelpers::fillHisto(effPlotCache.n_reco_times, RecoTimes);
  }
  
  void
  PhysicsPerformancePlotTool::fill(PhysicsPerformancePlotCache& effPlotCache,
				   const ActsFatras::Particle& truthParticle) const
  {
    const auto t_pT = truthParticle.transverseMomentum();
    const auto t_eta = eta(truthParticle.unitDirection());
    const auto t_phi = phi(truthParticle.unitDirection());

    this->fill(effPlotCache, t_pT, t_eta, t_phi);
  }

}
