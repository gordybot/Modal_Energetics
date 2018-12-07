ftn -o vmodes_feb.x vmodes_feb.f -traceback -g -O3 -fp-model source -warn nogeneral -convert big_endian -assume byterecl -mcmodel=small -fPIC

