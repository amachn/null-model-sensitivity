#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <vector>

using namespace Rcpp;

struct GcIdx {
  const Rbyte* words;
  const int* prefix;
  const int* ns;
  const int* ne;
  const int* nc;
  int nr;
  int len;
};

static std::vector<GcIdx> gc_unpack(const List& idx) {
  std::vector<GcIdx> result(idx.size());
  
  for (R_xlen_t i = 0; i < idx.size(); ++i) {
    List e = idx[i];
    RawVector w = e["words"]; IntegerVector pf = e["prefix"];
    IntegerVector ns = e["nstart"], ne = e["nend"], nc = e["ncum"];
    
    GcIdx g;
    g.words = w.begin(); g.prefix = pf.begin(); 
    g.ns = ns.begin(); g.ne = ne.begin(); g.nc = nc.begin(); g.nr = ns.size();
    g.len = as<int>(e["len"]);
    
    result[i] = g;
  }
  
  return result;
}

// return the number of G/C bases before position x
static inline int gc_before(const GcIdx& g, int x) {
  int word = x >> 6;
  int rem = x & 63;
  int count = g.prefix[word];
  
  if (rem) {
    uint64_t val;
    std::memcpy(&val, g.words + (size_t)word * 8, 8);
    count += __builtin_popcountll(val & ((1ULL << rem) - 1ULL));
  }
  
  return count;
}

// return the number of N bases before position x
static inline int n_before(const GcIdx& g, int x) {
  if (g.nr == 0) return 0;
  int j = (int)(std::lower_bound(g.ns, g.ns + g.nr, x) - g.ns);
  if (j == 0) return 0;
  int k = j - 1;
  return g.nc[k] + std::min(x, g.ne[k]) - g.ns[k];
}

// build one reusable chromosome index
// [[Rcpp::export]]
List gc_build_index(RawVector seq) {
  const int length = seq.size();
  const int nw = (length + 63) / 64;
  
  bool is_gc[256], is_n[256];
  for (int i = 0; i < 256; ++i) { is_gc[i] = false; is_n[i] = false; }
  
  is_gc[71] = is_gc[103] = is_gc[67] = is_gc[99] = true; // G, g, C, c
  is_n[78] = is_n[110] = true; // N, n
  
  std::vector<uint64_t> buf(nw, 0ULL);
  std::vector<int> ns, ne;
  
  int run_start = -1;
  const Rbyte* bases = seq.begin();
  
  for (int i = 0; i < length; ++i) {
    unsigned char base = bases[i];
    if (is_gc[base]) buf[i >> 6] |= (1ULL << (i & 63));
    
    if (is_n[base]) { if (run_start < 0) run_start = i; }
    else if (run_start >= 0) { 
      ns.push_back(run_start); ne.push_back(i); run_start = -1; 
    }
  }
  
  if (run_start >= 0) { ns.push_back(run_start); ne.push_back(length); }
  
  IntegerVector prefix(nw + 1);
  int acc = 0;
  
  for (int w = 0; w < nw; ++w) { 
    prefix[w] = acc; acc += __builtin_popcountll(buf[w]); 
  }
  
  prefix[nw] = acc;
  RawVector words((R_xlen_t)nw * 8);
  
  if (nw > 0) std::memcpy(words.begin(), buf.data(), (size_t)nw * 8);
  
  IntegerVector nstart(ns.begin(), ns.end()), nend(ne.begin(), ne.end());
  IntegerVector ncum(ns.size());
  int cum = 0;
  
  for (size_t k = 0; k < ns.size(); ++k) { ncum[k] = cum; cum += ne[k] - ns[k]; }
  
  return List::create(
    _["words"] = words,
    _["prefix"] = prefix,
    _["nstart"] = nstart,
    _["nend"] = nend,
    _["ncum"] = ncum,
    _["len"] = length
  );
}

// count G/C and N bases for 1-based inclusive intervals
// [[Rcpp::export]]
List gc_window_stats(List idx, IntegerVector chr_id, 
                     IntegerVector left, IntegerVector right) {
  std::vector<GcIdx> chromosomes = gc_unpack(idx);
  
  const R_xlen_t n = chr_id.size();
  NumericVector gc(n), nn(n);
  
  for (R_xlen_t i = 0; i < n; ++i) {
    const GcIdx& chromosome = chromosomes[chr_id[i] - 1];
    int start = left[i] - 1, end = right[i];
    
    gc[i] = gc_before(chromosome, end) - gc_before(chromosome, start);
    nn[i] = n_before(chromosome, end) - n_before(chromosome, start);
  }
  
  return List::create(_["gc"] = gc, _["n"] = nn);
}