#include "tbox/Allocator.h"

#include <cmath>
#include <deque>
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

constexpr std::size_t s_large_threshold = 4 * 1024 * 1024;
constexpr std::size_t s_n_large_blocks = 128;

std::size_t get_block_id(const std::size_t block_size) {
  return (block_size < 2 ? 0 : std::ilogb(block_size - 1) + 1);
}

/**
 * Manually manage the pools of smaller and larger objects.
 *
 * These objects must be deleted (and their stored buffers freed) in
 * Allocator::~Allocator(). At the same time, it is preferrable to not make them
 * class members so they can be used in static functions (and also so we do not
 * need to make them visible outside this translation unit), so they need to be
 * created, but not destroyed by static functions: hence, create them with
 * pointers and destroy them in Allocator::~Allocator().
 *
 * Rather than allocating on the heap, create these objects in static buffers to
 * decrease the chance a second exception (std::bad_alloc) is thrown if we
 * allocate during stack unwinding. std::deque's default ctor allocates memory
 * so we cannot work around that case, but since that's only for large buffers
 * it is unlikely to happen (e.g., we would have to run out of memory, throw an
 * exception, and then try to allocate a large array for the first time) during
 * error processing so, if that occurs, let the application terminate.
 */
std::vector<std::vector<void *>> *&
get_block_stacks()
{
  using vector_type = std::vector<std::vector<void *>>;
  alignas(vector_type) static char buffer[sizeof(vector_type)];
  static auto *s_block_stacks = new(buffer) vector_type();
  return s_block_stacks;
}

std::deque<std::pair<void*, std::size_t>> *&
get_large_blocks()
{
  using deque_type = std::deque<std::pair<void*, std::size_t>>;
  alignas(deque_type) static char buffer[sizeof(deque_type)];
  static auto *s_large_blocks = new(buffer) deque_type();
  return s_large_blocks;
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

  for (auto &pair : *get_large_blocks()) {
    std::free(pair.first);
  }

  // 1 of 2: we may run ~Allocator() before every ~Array() is run. To
  // avoid problems, clear data and set the boolean to false
  (*get_block_stacks()).~vector();
  get_block_stacks() = nullptr;
  (*get_large_blocks()).~deque();
  get_large_blocks() = nullptr;
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
  const std::size_t allocation_size = std::size_t(1) << block_id;
  TBOX_ASSERT(allocation_size >= n_bytes);
  if (allocation_size < s_large_threshold) {
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
  } else {
    auto &large_blocks = *get_large_blocks();
    auto it = large_blocks.end();
    for (std::size_t i = 0; i < large_blocks.size(); ++i) {
        if (large_blocks[i].second == n_bytes) {
            it = large_blocks.begin() + i;
            break;
        }
    }

    if (it == large_blocks.end()) {
        auto *block = std::malloc(n_bytes);
        return std::make_pair(true, block);
    } else {
        auto *block = it->first;
        large_blocks.erase(it);
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
    const auto block_id = get_block_id(n_bytes);
    const std::size_t allocation_size = std::size_t(1) << block_id;
    TBOX_ASSERT(allocation_size >= n_bytes);
    if (allocation_size < s_large_threshold) {
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
    } else {
      auto &large_blocks = *get_large_blocks();
      if (large_blocks.size() == s_n_large_blocks) {
        std::free(large_blocks.front().first);
        large_blocks.pop_front();
      }

      // same
      static std::size_t counter = 0;
      ++counter;
      if (counter % 1024 == 0) {
          const auto n_frees = large_blocks.size() / 4;
          for (std::size_t i = 0; i < n_frees; ++i) {
              std::free(large_blocks.front().first);
              large_blocks.pop_front();
          }
      }

      large_blocks.emplace_back(buffer, n_bytes);
    }
  }
}

}
}
