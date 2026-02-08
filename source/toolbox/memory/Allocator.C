#include <tbox/Allocator.h>

namespace SAMRAI {
namespace tbox {

bool Allocator::s_is_available = false;

std::vector<std::vector<void *>> Allocator::s_block_stacks;

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



}
}
