#include "write/allele_counter.h"
#include "utils/logger_public.h"
#include "utils/utils.h"

namespace tiledb::vcf {

AlleleCounter::AlleleCounter() {
}

AlleleCounter::~AlleleCounter() {
  if (dst_ != nullptr) {
    free(dst_);
  }
  finalize();
}

void AlleleCounter::create(
    Context& ctx, std::string root_uri, tiledb_filter_type_t checksum) {
  // Create filters
  Filter gzip(ctx, TILEDB_FILTER_GZIP);
  gzip.set_option(TILEDB_COMPRESSION_LEVEL, 9);

  FilterList filters_1(ctx);
  FilterList filters_2(ctx);
  filters_1.add_filter(gzip);
  filters_2.add_filter({ctx, TILEDB_FILTER_DOUBLE_DELTA});
  filters_2.add_filter(gzip);
  if (checksum) {
    filters_1.add_filter({ctx, checksum});
    filters_2.add_filter({ctx, checksum});
  }

  // Create schema and domain
  ArraySchema schema(ctx, TILEDB_SPARSE);
  Domain domain(ctx);
  schema.set_allows_dups(true);
  auto allele =
      Dimension::create(ctx, AC_ALLELE, TILEDB_STRING_ASCII, nullptr, nullptr);
  domain.add_dimension(allele);
  schema.set_domain(domain);
  schema.set_coords_filter_list(filters_1);
  schema.set_offsets_filter_list(filters_2);

  // Create attributes
  auto count = Attribute::create<int32_t>(ctx, AC_COUNT, filters_1);
  schema.add_attribute(count);

  // Create array
  // TODO: create utils function
  auto delim = utils::starts_with(root_uri, "tiledb://") ? '-' : '/';
  auto uri = utils::uri_join(root_uri, AC_URI, delim);
  Array::create(uri, schema);
}

void AlleleCounter::init(std::shared_ptr<Context> ctx, std::string root_uri) {
  ctx_ = ctx;
  auto uri = utils::uri_join(root_uri, AC_URI, '/');
  array_.reset(new Array(*ctx_, uri, TILEDB_WRITE));
  if (array_ == nullptr) {
    LOG_FATAL("AlleleCounter: error opening array '{}'", uri);
  }
}

void AlleleCounter::flush() {
  if (ac_count_.size() == 0) {
    return;
  }

  // Update results for the last locus before flushing
  update_results();

  LOG_DEBUG("AlleleCounter: flush {} records", ac_count_.size());
  LOG_DEBUG("AlleleCounter: allele = {}...", ac_allele_.substr(0, 64));

  Query query(*ctx_, *array_);

  auto st = query.set_layout(TILEDB_UNORDERED)
                .set_data_buffer(AC_ALLELE, ac_allele_)
                .set_offsets_buffer(AC_ALLELE, ac_allele_offsets_)
                .set_data_buffer(AC_COUNT, ac_count_)
                .submit();

  if (st != Query::Status::COMPLETE) {
    LOG_FATAL("AlleleCounter: error submitting TileDB write query");
  }

  ac_allele_.clear();
  ac_allele_offsets_.clear();
  ac_count_.clear();
}

void AlleleCounter::finalize() {
  flush();
  if (array_ != nullptr) {
    array_->close();
    array_ = nullptr;
  }
}

void AlleleCounter::update_results() {
  if (allele_count_.size() > 0) {
    for (auto& [allele, count] : allele_count_) {
      ac_allele_offsets_.push_back(ac_allele_.size());
      ac_allele_ += fmt::format("{}:{}", locus_, allele);
      ac_count_.push_back(count);
    }
    allele_count_.clear();
  }
}

void AlleleCounter::process(
    bcf_hdr_t* hdr,
    const std::string& sample_name,
    const std::string& contig,
    uint32_t pos,
    bcf1_t* rec) {
  // Check if locus has changed
  std::string locus = fmt::format("{}:{}", contig, pos);
  if (locus != locus_) {
    update_results();
    locus_ = locus;
  }

  // TODO: should we normalize REF, ALT alleles?
  int ngt = bcf_get_genotypes(hdr, rec, &dst_, &ndst_);
  for (int i = 0; i < ngt; i++) {
    // Skip missing and REF alleles
    if (bcf_gt_is_missing(dst_[i]) || bcf_gt_allele(dst_[i]) == 0) {
      continue;
    }
    std::string allele = rec->d.allele[bcf_gt_allele(dst_[i])];
    allele_count_[allele]++;
  }
}

}  // namespace tiledb::vcf
