#include <tbox/Allocator.h>

#include <cmath>
#include <deque>
#include <cstdlib>
#include <vector>

namespace SAMRAI {
namespace tbox {

namespace {
bool s_is_available = false;
std::vector<std::vector<void *>> s_block_stacks;
std::deque<std::pair<void*, std::size_t>> s_large_blocks;

constexpr std::size_t s_large_threshold = 4 * 1024 * 1024;
constexpr std::size_t s_n_large_blocks = 128;

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
  for (std::size_t i = 0; i < s_block_stacks.size(); ++i) {
     auto &block_stack = s_block_stacks[i];
     for (auto &block : block_stack) {
        std::free(block);
     }
  }

  for (auto &pair : s_large_blocks) {
     std::free(pair.first);
  }

  // 1 of 2: we may run ~Allocator() before every ~Array() is run. To
  // avoid problems, clear data and set the boolean to false
  s_block_stacks = {};
  s_large_blocks = {};
  s_is_available = false;
}

std::pair<bool, void *>
Allocator::internal_allocate(std::size_t n_bytes) {
  if (!s_is_available) {
     // In unusual circumstances (such as code running after main() finishes) we
     // may allocate memory after ~Allocator() is called: in that case, just
     // completely ignore the pool infrastructure
     return std::make_pair(true, std::malloc(n_bytes));
  }

  const auto block_id = get_block_id(n_bytes);
  const std::size_t allocation_size = 1 << block_id;
  TBOX_ASSERT(allocation_size >= n_bytes);
  if (allocation_size < s_large_threshold) {
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
  } else {
    auto it = s_large_blocks.end();
    for (std::size_t i = 0; i < s_large_blocks.size(); ++i) {
        if (s_large_blocks[i].second == n_bytes) {
            it = s_large_blocks.begin() + i;
            break;
        }
    }

    if (it == s_large_blocks.end()) {
        auto *block = std::malloc(n_bytes);
        return std::make_pair(true, block);
    } else {
        auto *block = it->first;
        s_large_blocks.erase(it);
        return std::make_pair(false, block);
    }
  }
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
