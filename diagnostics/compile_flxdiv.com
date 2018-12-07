ftn -o flx_diverg_221.x flx_diverg_221.f -traceback -g -O3 -fp-model source -warn nogeneral -convert big_endian -assume byterecl -mcmodel=small -fPIC

