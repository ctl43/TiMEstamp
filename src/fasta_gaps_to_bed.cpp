#include <Rcpp.h>
#include <fstream>
#include <string>
#include <vector>
#include <cctype>
#include <algorithm>

// Return basename without extension
static std::string basename_no_ext(const std::string& path) {
  size_t slash_pos = path.find_last_of("/\\");
  std::string base = (slash_pos == std::string::npos) ? path : path.substr(slash_pos + 1);
  size_t dot_pos = base.find_last_of('.');
  if (dot_pos != std::string::npos) {
    return base.substr(0, dot_pos);
  }
  return base;
}

// Trim trailing CR from lines (for CRLF files)
static inline void chomp_cr(std::string& s) {
  if (!s.empty() && s.back() == '\r') s.pop_back();
}

// [[Rcpp::export]]
void cxx_fasta_gaps_to_bed(Rcpp::CharacterVector fasta_files,
                           std::string bed_out,
                           int min_gap_inclusive = 10) {
  // sanitize threshold
  const int MIN_GAP = std::max(0, min_gap_inclusive);
  
  // Open output BED
  std::ofstream out(bed_out.c_str(), std::ios::out | std::ios::trunc);
  if (!out.is_open()) {
    Rcpp::stop("Cannot open output BED file: " + bed_out);
  }
  
  // Process each FASTA file
  for (int idx = 0; idx < fasta_files.size(); ++idx) {
    std::string fa_path = Rcpp::as<std::string>(fasta_files[idx]);
    std::ifstream in(fa_path.c_str());
    if (!in.is_open()) {
      out.close();
      Rcpp::stop("Cannot open FASTA file: " + fa_path);
    }
    
    const std::string file_tag = basename_no_ext(fa_path);
    
    std::string line;
    std::string chrom = "";        // current record name (from '>' header; first token)
    std::size_t pos = 0;           // 0-based position within current record
    bool in_gap = false;
    std::size_t gap_start = 0;
    
    auto emit_gap_if_long_enough = [&](std::size_t start, std::size_t end){
      std::size_t len = (end > start) ? (end - start) : 0;
      if (static_cast<int>(len) >= MIN_GAP) {
        // BED: chrom, start, end, name
        out << chrom << '\t' << start << '\t' << end << '\t' << file_tag << '\n';
      }
    };
    
    auto flush_gap = [&](){
      if (in_gap) {
        emit_gap_if_long_enough(gap_start, pos);
        in_gap = false;
      }
    };
    
    auto reset_record = [&](const std::string& new_chrom){
      // close any ongoing gap from previous record
      flush_gap();
      chrom = new_chrom;
      pos = 0;
      in_gap = false;
    };
    
    while (std::getline(in, line)) {
      chomp_cr(line);
      if (line.empty()) continue;
      
      if (line[0] == '>') {
        // Parse new header: take first token after '>'
        std::string hdr = line.substr(1);
        // trim leading spaces
        std::size_t i = 0;
        while (i < hdr.size() && std::isspace(static_cast<unsigned char>(hdr[i]))) ++i;
        std::size_t j = i;
        while (j < hdr.size() && !std::isspace(static_cast<unsigned char>(hdr[j]))) ++j;
        std::string new_chrom = hdr.substr(i, j - i);
        if (new_chrom.empty()) {
          in.close();
          out.close();
          Rcpp::stop("Malformed FASTA header (empty identifier) in file: " + fa_path);
        }
        reset_record(new_chrom);
        continue;
      }
      
      // Sequence line: iterate characters; count only non-whitespace
      for (char c : line) {
        if (std::isspace(static_cast<unsigned char>(c))) continue;
        
        if (c == '-') {
          if (!in_gap) {
            gap_start = pos;
            in_gap = true;
          }
          ++pos;
        } else {
          // non-gap base
          if (in_gap) {
            // close gap: [gap_start, pos)
            emit_gap_if_long_enough(gap_start, pos);
            in_gap = false;
          }
          ++pos;
        }
      }
    }
    
    // End of file: close any open gap and record
    flush_gap();
    in.close();
  }
  
  out.close();
}