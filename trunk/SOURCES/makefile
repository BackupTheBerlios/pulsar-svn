OBJECTS = common_module.o constantes_module.o ll.o pulsar.o subroutines.o

CC = gfortran

CFLAGS = -Wall 

pulsar.exe: $(OBJECTS)
	$(CC) -o pulsar.exe $(OBJECTS)

common_module.o: common_module.f95
	$(CC) $(CFLAGS) -c -o common_module.o common_module.f95
	
constantes_module.o: constantes_module.f95
	$(CC) $(CFLAGS) -c -o constantes_module.o constantes_module.f95

ll.o: ll.f95
	$(CC) $(CFLAGS) -c -o ll.o ll.f95

pulsar.o: pulsar.f95
	$(CC) $(CFLAGS) -c -o pulsar.o pulsar.f95

subroutines.o: subroutines.f95
	$(CC) $(CFLAGS) -c -o subroutines.o subroutines.f95

clean: 
	del *.o *.mod *.bak

.PHONY: clean