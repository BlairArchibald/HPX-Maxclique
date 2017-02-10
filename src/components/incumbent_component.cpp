#include "incumbent_component.hpp"

namespace globalBound
{
  void incumbent::updateBound(int size, std::set<int> members) {
    if (size > _size) {
      _size = size;
      _members = members;
    }
  }

  std::set<int> incumbent::getMembers() {
      return _members;
  }

  int incumbent::getBound() {
    return _size;
  }
}
HPX_REGISTER_COMPONENT_MODULE();

typedef hpx::components::component<globalBound::incumbent> incumbent_type;

HPX_REGISTER_COMPONENT(incumbent_type, incumbent);

HPX_REGISTER_ACTION(globalBound::incumbent::updateBound_action, incumbent_updateBound_action);
HPX_REGISTER_ACTION(globalBound::incumbent::getBound_action, incumbent_getBound_action);
HPX_REGISTER_ACTION(globalBound::incumbent::getMembers_action, incumbent_getMembers_action);
