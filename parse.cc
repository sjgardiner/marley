#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <regex>
#include <string>
#include <vector>
#include "ensdf_utils.hh"

class TEnsdfLevel;
class TEnsdfGamma;
class TEnsdfDecayScheme;

/// A class that represents an ENSDF nuclear energy level

class TEnsdfLevel {
  public:
    /// @param energy a string containing the
    /// level energy in keV
    /// @param jpi a string containing the spin
    /// and parity of the level (e.g., 0+)
    TEnsdfLevel(std::string energy = "0", std::string jpi = "");
    void add_gamma(const TEnsdfGamma& gamma);
    void clear_gammas();
    std::vector<TEnsdfGamma>* get_gammas();
    double get_numerical_energy() const;
    std::string get_string_energy() const;
    std::string get_spin_parity() const;
    void set_energy(std::string energy);
    void set_spin_parity(std::string jpi);

  private:
    std::string sEnergy;
    double fEnergy;
    std::string spin_parity;
    std::vector<TEnsdfGamma> gammas;
};

TEnsdfLevel::TEnsdfLevel(std::string energy, std::string jpi) {
  sEnergy = energy;
  fEnergy = std::stod(energy); // Converts the energy string into a double
  spin_parity = jpi;
}

void TEnsdfLevel::add_gamma(const TEnsdfGamma& gamma) {
  gammas.push_back(gamma);
}

double TEnsdfLevel::get_numerical_energy() const {
  return fEnergy; 
}

std::string TEnsdfLevel::get_string_energy() const {
  return sEnergy; 
}

void TEnsdfLevel::set_energy(std::string energy) {
  sEnergy = energy;
  fEnergy = std::stod(energy);
}

std::string TEnsdfLevel::get_spin_parity() const {
  return spin_parity;
}

void TEnsdfLevel::set_spin_parity(std::string jpi) {
  spin_parity = jpi;
}

void TEnsdfLevel::clear_gammas() {
  gammas.clear();
}

std::vector<TEnsdfGamma>* TEnsdfLevel::get_gammas() {
  return &gammas;
}


class TEnsdfGamma {
  public:
    TEnsdfGamma(double energy = 0, double ri = 0, TEnsdfLevel* start_level = nullptr);
    void set_start_level(TEnsdfLevel* start_level);
    void set_end_level(TEnsdfLevel* end_level);
    TEnsdfLevel* get_end_level() const;
    TEnsdfLevel* get_start_level() const;
    double get_energy() const;
    double get_ri() const;

  private:
    double fEnergy;
    double fRI;
    double fCC;
    double fTI; 
    TEnsdfLevel* pStartLevel;
    TEnsdfLevel* pEndLevel;
};

TEnsdfGamma::TEnsdfGamma(double energy, double ri, TEnsdfLevel* start_level) {
  fEnergy = energy;
  fRI = ri;
  fCC = 0;
  fTI = 0;
  pStartLevel = start_level;
  pEndLevel = nullptr;
}

void TEnsdfGamma::set_start_level(TEnsdfLevel* start_level) {
  pStartLevel = start_level;
}

void TEnsdfGamma::set_end_level(TEnsdfLevel* end_level) {
  pEndLevel = end_level;
}

TEnsdfLevel* TEnsdfGamma::get_start_level() const {
  return pStartLevel;
}

TEnsdfLevel* TEnsdfGamma::get_end_level() const {
  return pEndLevel;
}

double TEnsdfGamma::get_energy() const {
  return fEnergy;
}

double TEnsdfGamma::get_ri() const {
  return fRI;
}


class TEnsdfDecayScheme {
  public:
    TEnsdfDecayScheme(std::string nucid = "");
    std::string get_nuc_id() const;
    void set_nuc_id(std::string id);
    void add_level(const TEnsdfLevel level);
    std::map<std::string, TEnsdfLevel>* get_levels();
    TEnsdfLevel* get_level(std::string energy);
    std::vector<TEnsdfLevel*>* get_sorted_level_pointers();
    //void doCascade();

  private:
    std::string nuc_id;
    std::map<std::string, TEnsdfLevel> levels;
    std::vector<TEnsdfLevel*> pv_sorted_levels;
    static bool compare_level_energies(TEnsdfLevel* first,
      TEnsdfLevel* second);
};

TEnsdfDecayScheme::TEnsdfDecayScheme(std::string nucid) {
  nuc_id = nucid;
}

std::string TEnsdfDecayScheme::get_nuc_id() const {
  return nuc_id;
}

void TEnsdfDecayScheme::set_nuc_id(std::string id) {
  nuc_id = id;
}

bool TEnsdfDecayScheme::compare_level_energies(TEnsdfLevel* first,
  TEnsdfLevel* second)
{
  return first->get_numerical_energy() < second->get_numerical_energy();
}

void TEnsdfDecayScheme::add_level(TEnsdfLevel level) {
  // Add the level to the std::map of level objects. Use the
  // string version of its energy in keV as the key.
  std::string energy_string = level.get_string_energy();
  levels[energy_string] = level;

  // Get a pointer to the just-added level object
  TEnsdfLevel* p_level = &(levels[energy_string]);

  // Figure out where this level should go in the
  // vector of level pointers sorted by ascending energy
  std::vector<TEnsdfLevel*>::iterator
    insert_point = std::lower_bound(pv_sorted_levels.begin(),
    pv_sorted_levels.end(), p_level,
    compare_level_energies);

  // Add a pointer to the new level object at the
  // appropriate place in the energy-sorted vector
  pv_sorted_levels.insert(insert_point, p_level);
}


std::vector<TEnsdfLevel*>* TEnsdfDecayScheme::get_sorted_level_pointers() {
  return &pv_sorted_levels;
}


std::map<std::string, TEnsdfLevel>* TEnsdfDecayScheme::get_levels() {
  return &levels;
}

TEnsdfLevel* TEnsdfDecayScheme::get_level(std::string energy) {
  return &(levels.at(energy));
}


std::string process_continuation_records(std::ifstream &file_in,
  std::string &record, std::regex &rx_cont_record) {

  std::string line;
  while (!file_in.eof()) {
    // Get the next line of the file
    std::getline(file_in, line);
 
    // Check to see if the next line is a continuation record 
    if (std::regex_match(line, rx_cont_record)) {
      // If it is, add the next line to the current
      // record text.
      record += "\n" + line;
    }

    // If a non-continuation-record line is found,
    // stop reading from the file.
    else return line; 
  
  }

}

bool numeric_less_than(const std::string& s1, const std::string& s2) {
  return ensdf_utils::str_to_double(s1) < ensdf_utils::str_to_double(s2);
}





int main() {

  std::string nuc_id = " 40K ";
  std::string filename = "ensdf.040";

  // Create a decay scheme object to store data
  // imported from the ENSDF file
  TEnsdfDecayScheme decay_scheme(nuc_id);

  // Regular expressions for identifying ensdf record types
  //std::string generic_nuc_id = "^[[:alnum:] ]{5}";
  std::string primary_record = "[ 1]";
  std::string continuation_record = "[^ 1]";
  std::string record_data = ".{71}";
  
  std::regex rx_end_record("\\s*"); // Matches blank lines
  //std::regex rx_generic_primary_identification_record(nuc_id + primary_record + "   " + record_data);
  std::regex rx_primary_identification_record(nuc_id + primary_record
    + "   ADOPTED LEVELS, GAMMAS.{49}");
  std::regex rx_primary_level_record(nuc_id + primary_record + " L " + record_data);
  std::regex rx_continuation_level_record(nuc_id + continuation_record + " L " + record_data);
  std::regex rx_primary_gamma_record(nuc_id + primary_record + " G " + record_data);
  std::regex rx_continuation_gamma_record(nuc_id + continuation_record + " G " + record_data);

  // Open the ENSDF file for parsing
  std::ifstream file_in(filename);

  std::string line; // String to store the current line
                    // of the ENSDF file during parsing

  std::string record;

  bool found_decay_scheme = false; // Flag that indicates whether or not
                                   // gamma decay scheme data were found for the
                                   // current nuc_id.

  while (!file_in.eof()) {
    std::getline(file_in, line);
    if (std::regex_match(line, rx_primary_identification_record)) {
      found_decay_scheme = true;
      break;
    }
  }

  if (!found_decay_scheme) {
    std::cout << "Gamma decay scheme data (adopted levels, gammas) for " + nuc_id << std::endl;
    std::cout << "could not be found in the file " + filename << std::endl;
  }
  else {
    std::cout << "Gamma decay scheme data for " + nuc_id << " found. Using ENSDF dataset" << std::endl;
    std::cout << line << std::endl;
  }

  bool no_advance = false; // Flag that prevents advancing through
                           // the ENSDF file when a continuation record
                           // is found

  TEnsdfLevel* p_current_level = nullptr; // Pointer to the current level object
                                          // being filled with gamma ray data 

  while (!file_in.eof()) {
    // Get the next line of the file
    // unless we already did
    if (!no_advance) std::getline(file_in, line);
    no_advance = false;

    // Level Record
    if (std::regex_match(line, rx_primary_level_record)) {
      record = line;
      line = process_continuation_records(file_in, record, rx_continuation_level_record);
      no_advance = true;
        
      std::cout << "level record:" << std::endl << record << std::endl;

      // Extract the level energy (in keV) as a trimmed string from the ENSDF level record
      std::string level_energy = ensdf_utils::trim_copy(record.substr(9,10)); 

      // Also extract the spin-parity of the level
      std::string spin_parity = ensdf_utils::trim_copy(record.substr(21,18));

      // Add a new level object to our decay scheme using these data
      decay_scheme.add_level(TEnsdfLevel(level_energy, spin_parity));

      // Update the current level pointer
      p_current_level = decay_scheme.get_level(level_energy);
    }

    // Gamma Record
    else if (std::regex_match(line, rx_primary_gamma_record)) {
      record = line;
      line = process_continuation_records(file_in, record, rx_continuation_gamma_record);
      no_advance = true;

      std::cout << "gamma record:" << std::endl;

      // Extract the gamma ray's energy and relative (photon)
      // intensity from the ENSDF gamma record
      double gamma_energy = std::stod(record.substr(9,10));
      double gamma_ri = ensdf_utils::str_to_double(record.substr(21,8));

      // If this gamma belongs to a level record, then add its
      // data to the corresponding level object. Gammas that
      // have no assigned level appear before any level records,
      // so p_current_level will be a null pointer for them.
      if (p_current_level != nullptr) {
        p_current_level->add_gamma(TEnsdfGamma(gamma_energy, gamma_ri, p_current_level)); 
      }

      std::cout << record << std::endl;
      std::cout << "   energy: " << record.substr(9,10) << std::endl;
      std::cout << "   relative photon intensity: " << record.substr(21,7) << std::endl;
      std::cout << "   total conversion coefficient: " << record.substr(55,6) << std::endl;
      std::cout << "   relative total transition intensity: " << record.substr(64,9) << std::endl;
    }

    else if (std::regex_match(line, rx_end_record)) {
      std::cout << "Finished parsing gamma decay scheme data." << std::endl;
      break;
    }

  }

  // Now that we've loaded our decay scheme with all of the data,
  // we can fill in the end level pointers for each of our gamma
  // ray objects. This will enable the generator to descend through
  // the decay chains after only looking up the starting level.
  // Also check decay scheme object contents (for debugging purposes)
  std::cout << std::endl << std::endl << "Beginning decay scheme data check" << std::endl;
  std::vector<TEnsdfLevel*>* level_pointers = decay_scheme.get_sorted_level_pointers();
  // Use this vector of level energies to find the closest final level's index.
  // We'll push level energies onto it in ascending order to make the searching
  // as efficient as possible (we don't need to worry about negative-energy gamma
  // decays to higher levels)
  std::vector<double> level_energies;

  // Generate a report about the contents of the decay scheme object
  for(std::vector<TEnsdfLevel*>::iterator j = level_pointers->begin();
    j != level_pointers->end(); ++j)
  {
    
    std::cout << "Level at " << (*j)->get_string_energy() << std::endl;
    std::vector<TEnsdfGamma>* p_gammas = (*j)->get_gammas();
    double initial_level_energy = (*j)->get_numerical_energy();

    for(std::vector<TEnsdfGamma>::iterator k = p_gammas->begin();
      k != p_gammas->end(); ++k) 
    {
      double gamma_energy = k->get_energy();

      // Approximate the final level energy so we can search for the final level
      double final_level_energy = initial_level_energy - gamma_energy; 

      // Search for the corresponding energy using the vector of lower energies
      std::vector<double>::iterator p_final_level_energy = std::lower_bound(
        level_energies.begin(), level_energies.end(), final_level_energy); 

      // Determine the index of the final level energy appropriately
      int e_index = std::distance(level_energies.begin(), p_final_level_energy);
      if (e_index == level_energies.size()) {
        // The calculated final level energy is greater than
        // every energy in the vector. We will therefore assume
        // that the gamma decay takes us to the highest level in
        // the vector. Its index is given by one less than the
        // number of elements in the vector, so subtract one from
        // our previous result.
        --e_index;
      }
      else if (e_index > 0) {
        // If the calculated index does not correspond to the
        // first element, we still need to check which of the
        // two levels found (one on each side) is really the
        // closest. Do so and reassign the index if needed.
        if (std::abs(final_level_energy - level_energies[e_index])
          > std::abs(final_level_energy - level_energies[e_index - 1]))
        {
          --e_index;
        }
      }

      // Use the index to assign the appropriate end level pointer to this gamma
      k->set_end_level((*level_pointers)[e_index]);

      std::cout << "  has a gamma with energy " << k->get_energy();
      std::cout << " (transition to level at "
        << k->get_end_level()->get_string_energy() << " keV)" << std::endl;
      std::cout << "    and relative photon intensity " << k->get_ri() << std::endl; 
    }

    // Add the current level's energy to the level energies to use while
    // searching for gamma final levels
    level_energies.push_back(initial_level_energy);
  }

  std::cout << std::endl << "Completed decay scheme data check." << std::endl;

  file_in.close();

}
