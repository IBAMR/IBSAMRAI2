#ifndef included_tbox_ArrayAllocator
#define included_tbox_ArrayAllocator

#include "SAMRAI_config.h"

#include "tbox/Utilities.h"

#include <cmath>
#include <cstdlib>
#include <vector>
#include <type_traits>

namespace SAMRAI {
namespace tbox {
  /**
   * Internal allocation class. Responsible for managing a memory pool for
   * all Array instances.
   */
  class Allocator {
  public:
     Allocator();

     Allocator(const Allocator &) = delete;

     Allocator &operator=(const Allocator &) = delete;

    ~Allocator();

     static Allocator &getAllocator() {
      // we need a static instance so that we can deallocate arrays after main()
      // finishes
      static Allocator s_allocator;
      return s_allocator;
    }

  private:
    template <typename TYPE>
    static std::size_t &count()
    {
      static std::size_t s_count = 0;
      return s_count;
    }

    template <typename TYPE>
    static std::size_t get_block_id(const std::size_t block_size) {
      return (block_size < 2 ? 0 : std::ilogb(block_size * sizeof(TYPE) - 1) + 1);
    }

  public:
    template <typename TYPE>
    static std::size_t getNumberOfAllocations() {
       return count<TYPE>();
    }

    /**
     * Return a block of memory with enough space for @p block_size elements of
     * type TYPE. This block may either be newly allocated (via std::malloc())
     * or taken from a pool of previously allocated data.
     *
     * This function is not static to ensure it is only called after Allocator()
     * called.
     */
    template <typename TYPE>
    TYPE *allocate(const std::size_t block_size) {
      if (block_size == 0)
        return nullptr;

      const std::size_t block_id = get_block_id<TYPE>(block_size);
      const std::size_t allocation_size = 1 << block_id;
      TBOX_ASSERT(block_size <= allocation_size);
      // In unusual circumstances (such as code running after main() finishes)
      // we may allocate memory after ~Allocator() is called: in that case, just
      // completely ignore the pool infrastructure
      if (!s_is_available) {
         return static_cast<TYPE *>(std::malloc(allocation_size));
      }

      if (block_id >= s_block_stacks.size()) {
        s_block_stacks.resize(block_id + 1);
      }

      if (s_block_stacks[block_id].empty()) {
        auto block =
            static_cast<TYPE *>(std::malloc(allocation_size));
        s_block_stacks[block_id].reserve(s_block_stacks[block_id].capacity() +
                                         1);
        s_block_stacks[block_id].push_back(block);
        ++count<TYPE>();
      }

      auto *block = static_cast<TYPE *>(s_block_stacks[block_id].back());
      s_block_stacks[block_id].pop_back();
      if (!std::is_fundamental<TYPE>::value) {
        for (std::size_t k = 0; k < block_size; ++k) {
          new (&block[k]) TYPE;
        }
      }
      return block;
    }

    /**
     * Return an allocated block to the allocator.
     *
     * This function is static since it may be run after ~Allocator() is run.
     */
    template <typename TYPE>
    static void deallocate(TYPE *const block, const std::size_t block_size) {
      if (block_size == 0)
        return;

      if (!std::is_trivially_destructible<TYPE>::value) {
        for (std::size_t k = 0; k < block_size; ++k) {
          block[k].~TYPE();
        }
      }

      // 2 of 2: don't return memory to the Allocator if it has already been
      // destructed
      if (s_is_available) {
        TBOX_ASSERT(get_block_id<TYPE>(block_size) < s_block_stacks.size());
        s_block_stacks[get_block_id<TYPE>(block_size)].push_back(block);
      } else {
        std::free(block);
      }
    }

  private:
    static bool s_is_available;

    static std::vector<std::vector<void *>> s_block_stacks;
  };
}
}

#endif
