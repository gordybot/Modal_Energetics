ftn -o APE_may.x APE_may.f -traceback -g -O3 -fp-model source -warn nogeneral -convert big_endian -assume byterecl -mcmodel=small -fPIC

