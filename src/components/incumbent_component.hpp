#ifndef INCUMBENT_COMPONENT_HPP
#define INCUMBENT_COMPONENT_HPP

#include <hpx/hpx.hpp>
#include <hpx/include/components.hpp>
#include <hpx/include/serialization.hpp>

#include <utility>
#include <set>

namespace globalBound
{
  class incumbent : public hpx::components::locking_hook<
    hpx::components::component_base<incumbent>
    >
    {
    private:
      unsigned _size = 0;
      std::set<int> _members = {};
    public:
      void updateBound(int size, std::set<int> members);
      HPX_DEFINE_COMPONENT_ACTION(incumbent, updateBound);
      std::set<int> getMembers();
      HPX_DEFINE_COMPONENT_ACTION(incumbent, getMembers);
      int getBound();
      HPX_DEFINE_COMPONENT_ACTION(incumbent, getBound);
    };
}

HPX_REGISTER_ACTION_DECLARATION(globalBound::incumbent::updateBound_action, incumbent_updateBound_action);
HPX_REGISTER_ACTION_DECLARATION(globalBound::incumbent::getMembers_action, incumbent_getMembers_action);
HPX_REGISTER_ACTION_DECLARATION(globalBound::incumbent::getBound_action, incumbent_getBound_action);

#endif
