/**
 * @section LICENSE
 *
 * The MIT License
 *
 * @copyright Copyright (c) 2019 TileDB, Inc.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

#include "write/record_heap.h"

namespace tiledb {
namespace vcf {

void RecordHeap::clear() {
  while (!heap_.empty())
    heap_.pop();
}

bool RecordHeap::empty() const {
  return heap_.empty();
}

void RecordHeap::insert(
    VCF* vcf,
    NodeType type,
    bcf1_t* record,
    uint32_t sort_start_pos,
    uint32_t sample_id,
    bool end_node) {
  // Sanity check sort_start_pos >= start.
  if (sort_start_pos < (uint32_t)record->pos) {
    HtslibValueMem val;
    std::string contig(bcf_seqname(vcf->hdr(), record));
    std::string str_type = type == NodeType::Record ? "record" : "anchor";
    throw std::runtime_error(
        "Error inserting " + str_type + " '" + contig + ":" +
        std::to_string(record->pos + 1) + "-" +
        std::to_string(VCF::get_end_pos(vcf->hdr(), record, &val) + 1) +
        "' into ingestion heap from sample ID " + std::to_string(sample_id) +
        "; sort start position " + std::to_string(sort_start_pos + 1) +
        " cannot be less than start.");
  }

  auto node = std::unique_ptr<Node>(new Node);
  node->vcf = vcf;
  node->type = type;
  node->record = record;
  node->sort_start_pos = sort_start_pos;
  node->sample_id = sample_id;
  node->end_node = end_node;
  heap_.push(std::move(node));
}

const RecordHeap::Node& RecordHeap::top() const {
  return *heap_.top();
}

void RecordHeap::pop() {
  heap_.pop();
}

}  // namespace vcf
}  // namespace tiledb
