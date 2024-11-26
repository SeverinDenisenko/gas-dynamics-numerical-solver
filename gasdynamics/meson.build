gasdynamicslib_include = include_directories('include')
gasdynamicslib_deps = declare_dependency(include_directories : gasdynamicslib_include)

gasdynamicslib = static_library('gasdynamicslib', sources : ['solver.cpp'], include_directories : gasdynamicslib_include, dependencies : [gasdynamicslib_deps])
