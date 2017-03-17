default:
	make prep && FoBiS.py build -t src/hci.f90 -j 2

prep: src/linked_list_module.f90 src/sort_module.f90 src/hash_table_module.f90

src/linked_list_module.f90: src/linked_list_module.fpp
	fypp src/linked_list_module.fpp src/linked_list_module.f90

src/sort_module.f90: src/sort_module.fpp
	fypp src/sort_module.fpp src/sort_module.f90

src/hash_table_module.f90: src/hash_table_module.fpp
	fypp src/hash_table_module.fpp src/hash_table_module.f90

clean:
	FoBiS.py clean
