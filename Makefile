default:
	FoBiS.py build -t hci.f90

test:
	FoBiS.py build && rm -rf tests && mkdir tests && mv *_test tests/ && run-parts tests -v	
