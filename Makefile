default:
	make prep && FoBiS.py build -t src/hci.f90

prep: 
	fypp src/linked_list_module.fpp src/linked_list_module.f90
	fypp src/utilities_module.fpp src/utilities_module.f90

clean:
	FoBiS.py clean
