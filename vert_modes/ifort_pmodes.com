ifort -o pmodes_apr.x pmodes_apr.f -traceback -g -O3 -fp-model source -warn nogeneral -convert big_endian -assume byterecl -mcmodel=medium -fPIC
