default:
	FoBiS.py build -t hci.f90

pp:
	fypp src/linked_list_module.fpp src/linked_list_module.f90
	fypp src/utilities_module.fpp src/utilities_module.f90

test:
	FoBiS.py build && rm -rf tests && mkdir tests && mv *_test tests/ && run-parts tests -v	
