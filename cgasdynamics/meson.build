cgasdynamicslib_include = include_directories('include')
cgasdynamicslib_deps = declare_dependency(include_directories : cgasdynamicslib_include)

cgasdynamicslib = shared_library(
	'cgasdynamicslib',
	sources : ['solver.cpp'],
	link_with : [gasdynamicslib, riemann_solver_lib],
	dependencies : [cgasdynamicslib_deps, gasdynamicslib_deps, riemann_solver_lib_deps]
)
