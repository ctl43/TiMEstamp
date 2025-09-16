#include <Rcpp.h>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <string>
#include <filesystem>
#include <algorithm>

using namespace Rcpp;
using namespace std;
namespace fs = std::filesystem;

// Extract species name from "species.chromosome"
string extract_species(const string& species_line) {
  size_t dot_pos = species_line.find('.');
  if (dot_pos != string::npos) return species_line.substr(0, dot_pos);
  return species_line;
}

// Extract chromosome name from "species.chromosome"
string extract_chromosome(const string& species_line) {
  size_t dot_pos = species_line.find('.');
  if (dot_pos != string::npos) return species_line.substr(dot_pos + 1);
  return species_line; // if no dot, return as-is
}

// Remove columns where the reference has '-'
string remove_insertions(const string& reference, const string& sequence) {
  string cleaned;
  cleaned.reserve(sequence.size());
  for (size_t i = 0; i < reference.size(); ++i) {
    if (reference[i] != '-') cleaned += (i < sequence.size()) ? sequence[i] : '-';
  }
  return cleaned;
}

// Count non-gap characters in reference
size_t get_adjusted_length(const string& reference) {
  return count_if(reference.begin(), reference.end(), [](char c){ return c != '-'; });
}

// Approximate memory usage of buffers
size_t estimate_memory(const unordered_map<string, string>& buffer) {
  size_t total = 0;
  for (const auto& kv : buffer) total += kv.second.size();
  return total;
}

// Flush per-species buffers to disk
void flush_buffer(unordered_map<string, ofstream>& species_files,
                  unordered_map<string, string>& buffer) {
  for (auto& kv : buffer) species_files[kv.first] << kv.second;
  buffer.clear();
}

// Read chrom.sizes into a map<chrom, size>
unordered_map<string, long long> read_chrom_sizes(const string& chrom_size_file) {
  ifstream in(chrom_size_file);
  if (!in.is_open()) stop("Failed to open chrom.sizes file: " + chrom_size_file);
  unordered_map<string, long long> sizes;
  string chrom;
  long long size;
  while (in >> chrom >> size) {
    sizes[chrom] = size;
  }
  in.close();
  if (sizes.empty()) stop("chrom.sizes file appears empty.");
  return sizes;
}

// [[Rcpp::export]]
void cxx_convert_maf_to_fasta(std::string maf_file,
                                 std::string species_file,
                                 std::string chrom_size_file,
                                 std::string output_folder,
                                 std::string reference_species = "",
                                 size_t buffer_limit_mb = 100) {
  
  // Read species list
  ifstream species_in(species_file);
  if (!species_in.is_open()) stop("Failed to open species list file.");
  unordered_set<string> species_set;
  vector<string> species_list;
  string species_name;
  while (getline(species_in, species_name)) {
    if (!species_name.empty()) {
      species_set.insert(species_name);
      species_list.push_back(species_name);
    }
  }
  species_in.close();
  if (species_list.empty()) stop("Species list is empty.");
  
  // Read chrom.sizes
  auto chrom_sizes = read_chrom_sizes(chrom_size_file);
  
  // Ensure output directory
  if (!fs::exists(output_folder)) fs::create_directories(output_folder);
  
  // Open per-species FASTA (overwrite) — headers will be written later
  unordered_map<string, ofstream> species_files;
  for (const auto& sp : species_list) {
    string filename = output_folder + "/" + sp + ".fa";
    species_files[sp].open(filename, ios::out);
    if (!species_files[sp].is_open()) stop("Failed to create output file for species: " + sp);
  }
  
  // Buffers
  unordered_map<string, string> alignment_buffer;
  size_t buffer_limit_bytes = buffer_limit_mb * 1024 * 1024;
  
  // Open MAF
  ifstream maf_in(maf_file);
  if (!maf_in.is_open()) stop("Failed to open MAF file.");
  
  // State for current block
  string line;
  unordered_map<string, string> temp_alignments;
  string temp_reference_sequence;
  
  // Track current ungapped reference offset since start of chromosome
  long long ref_ungapped_offset = 0;
  // Size for the reference of the current block
  long long current_block_ref_size  = -1;
  
  // Reference chrom info
  string reference_chrom = "";
  long long reference_chrom_size = -1;
  bool reference_chrom_locked = false;
  bool header_written = false;
  
  auto write_headers_if_needed = [&](){
    if (!header_written && !reference_chrom.empty()) {
      for (auto& kv : species_files) {
        kv.second << ">" << reference_chrom << "\n";
      }
      header_written = true;
    }
  };
  
  // Insert inter-block gap: reference gets 'N', non-reference gets '-'
  auto add_interblock_gap = [&](long long gap){
    if (gap <= 0) return;
    string pad_ref(static_cast<size_t>(gap), 'N');
    string pad_nonref(static_cast<size_t>(gap), '-');
    for (const auto& sp : species_list) {
      alignment_buffer[sp] += (sp == reference_species ? pad_ref : pad_nonref);
    }
    // Ensure headers exist before any possible flush
    write_headers_if_needed();
    if (estimate_memory(alignment_buffer) >= buffer_limit_bytes) {
      flush_buffer(species_files, alignment_buffer);
    }
    ref_ungapped_offset += gap;
  };
  
  auto finalize_block = [&](){
    if (temp_reference_sequence.empty()) return;
    
    // Adjusted length equals number of non-gap characters in reference
    size_t adjusted_len = get_adjusted_length(temp_reference_sequence);
    
    // Append sequences with ref-insertions removed; fill missing with '-' of adjusted_len
    for (const auto& sp : species_list) {
      if (temp_alignments.find(sp) != temp_alignments.end()) {
        alignment_buffer[sp] += remove_insertions(temp_reference_sequence, temp_alignments[sp]);
      } else {
        alignment_buffer[sp] += string(adjusted_len, '-');
      }
    }
    
    // Ensure headers exist before any possible flush
    write_headers_if_needed();
    if (estimate_memory(alignment_buffer) >= buffer_limit_bytes) {
      flush_buffer(species_files, alignment_buffer);
    }
    
    // Advance ungapped reference offset by the block's ungapped size
    if (current_block_ref_size >= 0) {
      ref_ungapped_offset += current_block_ref_size;
    } else {
      ref_ungapped_offset += static_cast<long long>(adjusted_len);
    }
    
    // Clear block state
    temp_alignments.clear();
    temp_reference_sequence.clear();
    current_block_ref_size  = -1;
  };
  
  // Parse MAF
  while (getline(maf_in, line)) {
    if (line.empty() || line[0] == '#') continue;
    
    if (line[0] == 'a') {
      // Starting a new block: finalize previous one
      finalize_block();
      continue;
    }
    
    if (line[0] == 's') {
      istringstream iss(line);
      string type, species_line, sequence;
      long long start;
      long long size;
      char strand;
      long long srcSize;
      iss >> type >> species_line >> start >> size >> strand >> srcSize >> sequence;
      
      string sp = extract_species(species_line);
      
      // Determine reference species if not provided (first 's' seen)
      if (reference_species.empty()) reference_species = sp;
      
      // If this is the first time we see the reference 's' line, lock reference chromosome
      if (sp == reference_species && temp_reference_sequence.empty() && !reference_chrom_locked) {
        reference_chrom = extract_chromosome(species_line);
        auto it = chrom_sizes.find(reference_chrom);
        if (it == chrom_sizes.end()) {
          stop("Reference chromosome " + reference_chrom + " not found in chrom.sizes.");
        }
        reference_chrom_size = it->second;
        reference_chrom_locked = true;
        
        // Now we know the header: write it once to all files
        write_headers_if_needed();
      }
      
      // If this line is the reference for the current block, handle inter-block gap first
      if (sp == reference_species && temp_reference_sequence.empty()) {
        // Leading gap up to 'start' (ungapped ref coordinates)
        add_interblock_gap(start - ref_ungapped_offset);
        // Record block ref size (ungapped)
        current_block_ref_size  = size;
      }
      
      // Store the raw aligned sequence for this species for this block
      temp_alignments[sp] = sequence;
      
      // Capture the reference aligned sequence once per block
      if (sp == reference_species && temp_reference_sequence.empty()) {
        temp_reference_sequence = sequence;
      }
    }
    
    // Ignore other line types (e.g., 'i', 'e', 'q')
  }
  
  // Final block at EOF
  finalize_block();
  
  // Tail padding to the end of the reference chromosome per chrom.sizes
  if (!reference_chrom_locked) {
    // No reference encountered: cannot set header or pad
    for (auto& kv : species_files) kv.second.close();
    maf_in.close();
    stop("No reference block encountered. Cannot determine reference chromosome for FASTA header and tail padding.");
  }
  long long tail_gap = reference_chrom_size - ref_ungapped_offset;
  if (tail_gap > 0) {
    add_interblock_gap(tail_gap); // same rule: N for ref, '-' for others
  }
  
  // One more time: ensure headers exist before final flush (should already be true)
  write_headers_if_needed();
  
  // Final flush
  if (!alignment_buffer.empty()) flush_buffer(species_files, alignment_buffer);
  
  // Close files
  for (auto& kv : species_files) kv.second.close();
  maf_in.close();
  
  Rcout << "MAF processing completed. FASTA headers set to " << reference_chrom
        << ". Inter-block gaps padded (ref=N, others='-'). "
        << "All sequences extended to chrom length " << reference_chrom_size
        << ". Saved to: " << output_folder << "\n";
}


