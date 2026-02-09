#include <tbox/Allocator.h>

#include <cmath>
#include <cstdlib>
#include <vector>

namespace SAMRAI {
namespace tbox {

namespace {
std::vector<std::vector<void *>> s_block_stacks;

bool s_is_available = false;

std::size_t get_block_id(const std::size_t block_size) {
  return (block_size < 2 ? 0 : std::ilogb(block_size - 1) + 1);
}
}

Allocator &
Allocator::getAllocator() {
  // we need a static instance so that we can deallocate arrays after main()
  // finishes
  static Allocator s_allocator;
  return s_allocator;
}

Allocator::Allocator() {
   s_is_available = true;
}

Allocator::~Allocator() {
  for (auto &block_stack : s_block_stacks) {
     for (auto &block : block_stack) {
        std::free(block);
     }
  }
  // 1 of 2: we may run ~Allocator() before every ~Array() is run. To
  // avoid problems, clear data and set the boolean to false
  s_block_stacks = {};
  s_is_available = false;
}

std::pair<bool, void *>
Allocator::internal_allocate(std::size_t n_bytes) {
  const auto block_id = get_block_id(n_bytes);
  const std::size_t allocation_size = 1 << block_id;
  TBOX_ASSERT(allocation_size >= n_bytes);
  // In unusual circumstances (such as code running after main() finishes) we
  // may allocate memory after ~Allocator() is called: in that case, just
  // completely ignore the pool infrastructure
  if (!s_is_available) {
    return std::make_pair(true, std::malloc(allocation_size));
  }

  if (block_id >= s_block_stacks.size()) {
    s_block_stacks.resize(block_id + 1);
  }

  bool new_allocation = false;
  if (s_block_stacks[block_id].empty()) {
    auto block = std::malloc(allocation_size);
    TBOX_ASSERT(block != nullptr);
    new_allocation = true;
    s_block_stacks[block_id].push_back(block);
  }

  auto *block = s_block_stacks[block_id].back();
  s_block_stacks[block_id].pop_back();
  return std::make_pair(new_allocation, block);
}

void
Allocator::internal_deallocate(void *buffer, std::size_t n_bytes) {
  if (!s_is_available) {
    // 2 of 2: don't return memory to the Allocator if it has already been
    // destructed
    std::free(buffer);
  } else {
    s_block_stacks[get_block_id(n_bytes)].push_back(buffer);
  }
}

}
}
