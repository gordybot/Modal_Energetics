ftn -o extract_rho_221.x extract_rho_221.f -traceback -g -O3 -fp-model source -warn nogeneral -convert big_endian -assume byterecl -mcmodel=small -fPIC

