project('gasdynamics', 'cpp', 'c', default_options : ['warning_level=3', 'cpp_std=c++23'])

riemann_solver = subproject('riemann-problem-solver')
riemann_solver_lib = riemann_solver.get_variable('riemann_solver_lib')
riemann_solver_lib_deps = riemann_solver.get_variable('riemann_solver_lib_deps')

subdir('gasdynamics')
subdir('cgasdynamics')

fs = import('fs')
python_target = meson.source_root() + '/pygasdynamics/' + fs.name(cgasdynamicslib.full_path())
run_target(
	'python',
	command : ['cp', cgasdynamicslib.path(), python_target]
)
