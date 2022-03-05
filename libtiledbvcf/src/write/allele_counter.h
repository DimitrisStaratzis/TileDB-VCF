#ifndef TILEDB_VCF_ALLELE_COUNTER_H
#define TILEDB_VCF_ALLELE_COUNTER_H

#include <map>
#include <string>
#include <vector>

#include <htslib/vcf.h>
#include <tiledb/tiledb>

namespace tiledb::vcf {

class IngestionTask {};

class AlleleCounter : public IngestionTask {
 public:
  AlleleCounter();

  ~AlleleCounter();

  // create array
  static void create(
      Context& ctx, std::string root_uri, tiledb_filter_type_t checksum);

  // open array
  void init(std::string root_uri);

  void process(
      bcf_hdr_t* hdr,
      const std::string& sample_name,
      const std::string& contig,
      uint32_t pos,
      bcf1_t* record);

  // create query, write, finalize
  void flush();

  void finalize();

 private:
  inline static const std::string AC_ALLELE = "allele";
  inline static const std::string AC_COUNT = "count";
  inline static const std::string AC_URI = "allele_count";

  std::unique_ptr<Context> ctx_ = nullptr;
  std::unique_ptr<Array> array_ = nullptr;
  std::map<std::string, int> allele_count_;
  std::string locus_;

  std::string ac_allele_;
  std::vector<uint64_t> ac_allele_offsets_;
  std::vector<uint32_t> ac_count_;

  // reusable htslib buffer for bcf_get_* functions
  int* dst_ = nullptr;

  // reusable htslib buffer size for bcf_get_* functions
  int ndst_ = 0;
};

}  // namespace tiledb::vcf

#endif
