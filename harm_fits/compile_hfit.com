ftn -o charmfit_uiso.x charmfit_uiso.f -traceback -g -O3 -fp-model source -warn nogeneral -convert big_endian -assume byterecl -mcmodel=small -fPIC
