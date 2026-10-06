#ifndef included_tbox_Allocator
#define included_tbox_Allocator

#include "SAMRAI_config.h"

#include "tbox/Utilities.h"

#include <cstdlib>
#include <new>
#include <utility>
#include <type_traits>

namespace SAMRAI {
namespace tbox {
  /**
   * Internal allocation class. Responsible for managing a memory pool for
   * all Array instances.
   */
  class Allocator {
  public:
     static Allocator &getAllocator();

  private:
     /**
      * Constructor. As there is only one Allocator at a time this is private.
      */
     Allocator();

     Allocator(const Allocator &) = delete;

     Allocator &operator=(const Allocator &) = delete;

    ~Allocator();

    template <typename TYPE>
    static std::size_t &count()
    {
      static std::size_t s_count = 0;
      return s_count;
    }

    /**
     * Internal allocation function.
     *
     * @return A pair whose first value indicates whether or not the memory was
     * allocated (i.e., `true` for a new allocation, `false` for reusing an
     * existing allocation), and whose second value is a pointer to at least @p
     * n_bytes of storage.
     */
    static
    std::pair<bool, void *>
    internal_allocate(std::size_t n_bytes);

    /**
     * Internal deallocation function.
     *
     * The second size must be the same as the value provided to
     * internal_allocate().
     */
    static
    void
    internal_deallocate(void *buffer, std::size_t n_bytes);

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

      auto pair = internal_allocate(block_size * sizeof(TYPE));
      if (pair.first) {
        ++count<TYPE>();
      }
      auto* block = static_cast<TYPE *>(pair.second);
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

      internal_deallocate(block, block_size * sizeof(TYPE));
    }
  };
}
}

#endif
