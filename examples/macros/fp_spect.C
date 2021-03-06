#include <iostream>
#include <string>

void fp_spect(const std::string& file_name, int pdg) {

  marley::MacroEventFileReader reader( file_name );
  marley::Event ev;

  long num_events = 0;

  std::vector<double> KE_vec;

  while ( reader >> ev ) {

    if (num_events % 1000 == 0) std::cout << "Event " << num_events << '\n';

    size_t num_finals = ev.final_particle_count();

    for (size_t f = 0u; f < num_finals; ++f) {
      const marley::Particle& fp = ev.final_particle( f );
      if ( fp.pdg_code() == pdg ) {
        KE_vec.push_back( fp.kinetic_energy() );
      }
    }

    ++num_events;
  }

  double KE_max = -1e30;
  double KE_min = 1e30;
  for (size_t k = 0; k < KE_vec.size(); ++k) {
    double ke = KE_vec.at(k);
    if (ke > KE_max) KE_max = ke;
    else if (ke < KE_min) KE_min = ke;
  }

  TString title_str;
  title_str.Form("kinetic energies for pdg = %d; kinetic energy (MeV);"
    "events", pdg);

  TH1D* KEs = new TH1D("KEs", title_str.Data(), 100, KE_max, KE_min);

  // Prevent this histogram from being automatically deleted when the
  // current TFile is closed
  KEs->SetDirectory(NULL);

  for (size_t j = 0u; j < KE_vec.size(); ++j) {
    KEs->Fill( KE_vec.at(j) );
  }

  TCanvas* c = new TCanvas;
  c->cd();

  gStyle->SetOptStat();

  KEs->SetStats(false);
  KEs->SetLineColor(kBlue);
  KEs->SetLineWidth(2);
  KEs->Draw("hist");

  std::cout << "Found " << KE_vec.size() << " particles with"
    << " pdg = " << pdg << " in " << num_events << " events\n";
}
