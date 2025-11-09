# Macro definitions
FC = ifort
LDFLAGS = -O4 -mkl
FFLAGS = -qopenmp -extend-source
# ==============================
PES=bkmp2.f

SRC = machina_basic.f90 globalPara.f90 ${PES} potent.f90 function.f basis.f90 initWP.f90 main.f90
OBJ = $(SRC:.f90=.o)
OBJ := $(OBJ:.f=.o)
TARGET = scatt

$(TARGET): $(OBJ)
	$(FC) $(OBJ) -o $@ $(LDFLAGS)

%.o: %.f90
	$(FC) $(FFLAGS) -c $<

%.o: %.f
	$(FC) $(FFLAGS) -c $<

clean:
	rm -f *.o *.mod $(TARGET)
# ==============================