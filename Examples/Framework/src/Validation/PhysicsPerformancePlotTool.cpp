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
  
}
