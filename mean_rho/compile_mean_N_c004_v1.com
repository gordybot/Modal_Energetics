ftn -o mean_N_c004_v1.x mean_N_c004_v1.f  -traceback -g -O3 -fp-model source -warn nogeneral -convert big_endian -assume byterecl -mcmodel=small -fPIC
