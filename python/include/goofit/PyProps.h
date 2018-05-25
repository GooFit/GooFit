
#define ADD_PROP(name, get, set, parent)                                                                               \
    .def_property(#name, &parent::get, &parent::set).def(#get, &parent::get).def(#set, &parent::set)

#define ADD_PROP_RO(name, get, parent) .def_property(#name, &parent::get, nullptr).def(#get, &parent::get)

#define ADD_PROP_WO(name, set, parent) .def_property(#name, nullptr, &parent::set).def(#set, &parent::set)
