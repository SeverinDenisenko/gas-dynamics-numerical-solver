project('gasdynamics', 'cpp', 'c', default_options : ['warning_level=3', 'cpp_std=c++23'])

subdir('gasdynamics')
subdir('cgasdynamics')

fs = import('fs')
python_target = meson.source_root() + '/pygasdynamics/' + fs.name(cgasdynamicslib.full_path())
run_target('python', command : ['cp', cgasdynamicslib.path(), python_target])
