default:
	make prep && FoBiS.py build -t src/hci.f90 -j 2

prep: src/linked_list_module.f90 src/utilities_module.f90

src/linked_list_module.f90: src/linked_list_module.fpp
	fypp src/linked_list_module.fpp src/linked_list_module.f90

src/utilities_module.f90: src/utilities_module.fpp
	fypp src/utilities_module.fpp src/utilities_module.f90

clean:
	FoBiS.py clean
