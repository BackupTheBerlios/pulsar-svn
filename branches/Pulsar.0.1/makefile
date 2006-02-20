OBJECTS = common_module.o constantes_module.o ll.o pulsar.o subroutines.o

F95 = gfortran

CFLAGS = -Wall -O2

pulsar.exe: $(OBJECTS)
	$(F95) -o pulsar.exe $(OBJECTS)

common_module.o: common_module.f95
	$(F95) $(CFLAGS) -c -o common_module.o common_module.f95
	
constantes_module.o: constantes_module.f95
	$(F95) $(CFLAGS) -c -o constantes_module.o constantes_module.f95

ll.o: ll.f95
	$(F95) $(CFLAGS) -c -o ll.o ll.f95

pulsar.o: pulsar.f95
	$(F95) $(CFLAGS) -c -o pulsar.o pulsar.f95

subroutines.o: subroutines.f95
	$(F95) $(CFLAGS) -c -o subroutines.o subroutines.f95

clean: 
	del pulsar.exe *.o *.mod *.bak

.PHONY: clean