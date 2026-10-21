#include "tbox/Allocator.h"
#include "tbox/Utilities.h"

#include <cmath>
#include <cstdlib>
#include <vector>

namespace SAMRAI {
namespace tbox {

namespace {
/**
 * Whether or not the Allocator is available.
 *
 * In unusual circumstances (e.g., tbox::Array objects destroyed after main()
 * finishes) we may need to deallocate memory after the Allocator has been
 * destroyed. In that case simply deallocate the memory.
 */
static bool s_is_available = false;

std::size_t get_block_id(const std::size_t block_size) {
  return (block_size < 2 ? 0 : std::ilogb(block_size - 1) + 1);
}

/**
 * Manually manage the pools of blocks, one for each power of two.
 *
 * The pools must be deleted (and their stored buffers freed) in
 * Allocator::~Allocator(). At the same time, it is preferrable to not make them
 * class members so they can be used in static functions (and also so we do not
 * need to make them visible outside this translation unit), so they need to be
 * created, but not destroyed by static functions: hence, create them with
 * pointers and destroy them in Allocator::~Allocator().
 *
 * Rather than allocating on the heap, create them in a static buffer to
 * decrease the chance a second exception (std::bad_alloc) is thrown if we
 * allocate during stack unwinding.
 */
std::vector<std::vector<void *>> *&
get_block_stacks()
{
  using vector_type = std::vector<std::vector<void *>>;
  alignas(vector_type) static char buffer[sizeof(vector_type)];
  static auto *s_block_stacks = new(buffer) vector_type();
  return s_block_stacks;
}
}

Allocator &
Allocator::getAllocator() {
  // We create exactly one Allocator so that all allocations are pooled.
  static Allocator s_allocator;
  return s_allocator;
}

Allocator::Allocator() {
   s_is_available = true;
}

Allocator::~Allocator() {
  for (auto &block_stack : *get_block_stacks()) {
    for (auto &block : block_stack) {
      std::free(block);
     }
  }

  // 1 of 2: we may run ~Allocator() before every ~Array() is run. To
  // avoid problems, clear data and set the boolean to false
  (*get_block_stacks()).~vector();
  get_block_stacks() = nullptr;
  s_is_available = false;
}

std::pair<bool, void *>
Allocator::internal_allocate(std::size_t n_bytes) {
  if (!s_is_available) {
     // In unusual circumstances (such as code running after main() finishes) we
     // may allocate memory after ~Allocator() is called: in that case, just
     // completely ignore the pool infrastructure
     auto *block = std::malloc(n_bytes);
     TBOX_ASSERT(block != nullptr);
     return std::make_pair(true, block);
  }

  const auto block_id = get_block_id(n_bytes);
  const std::size_t allocation_size = std::size_t(1) << block_id;
  TBOX_ASSERT(allocation_size >= n_bytes);
  auto &block_stacks = *get_block_stacks();
  if (block_id >= block_stacks.size()) {
    block_stacks.resize(block_id + 1);
  }

  bool new_allocation = false;
  if (block_stacks[block_id].empty()) {
    auto block = std::malloc(allocation_size);
    TBOX_ASSERT(block != nullptr);
    new_allocation = true;
    block_stacks[block_id].push_back(block);
  }

  auto *block = block_stacks[block_id].back();
  block_stacks[block_id].pop_back();
  return std::make_pair(new_allocation, block);
}

void
Allocator::internal_deallocate(void *buffer, std::size_t n_bytes) {
  if (!s_is_available) {
    // 2 of 2: don't return memory to the Allocator if it has already been
    // destructed
    std::free(buffer);
  } else {
    const auto block_id = get_block_id(n_bytes);
    const std::size_t allocation_size = std::size_t(1) << block_id;
    TBOX_ASSERT(allocation_size >= n_bytes);
    auto &block_stacks = *get_block_stacks();
    block_stacks[get_block_id(n_bytes)].push_back(buffer);
    // Every so often, partially clear the cache to ensure massive allocations
    // aren't sticking around. This is based on some profiling which shows
    // that, over 100 time steps and 8 processors, we allocate about 500k
    // small things and 20k large things
    static std::size_t counter = 0;
    ++counter;
    if (counter % 65536 == 0) {
        for (auto &block_stack : block_stacks) {
            const auto n_frees = 3 * block_stack.size() / 4;
            for (std::size_t i = 0; i < n_frees; ++i) {
                std::free(block_stack[i]);
            }
            block_stack.erase(block_stack.begin(),
                              block_stack.begin() + n_frees);
        }
    }
  }
}

}
}
