# the following is the library made by compiler, and contains all the files.
# when you run make that will create lib.a

Exec = TaxSim

ObjDir = ./Objects
CodeDir = .

Objs = NPSOL.o CStoGrowthCom.o CStoGrowthOpt.o CStoGrowthDP.o
testObjs = NPSOL.o testNPSOL.o

Compiler = gfortran

# Fortran Compiler options
# next we choose the options for the chosen compiler

CompilerFlags = -O3 -fbounds-check -Wunused -Wunused-parameter

StaticOpt = -static-libgfortran
legacy = -std=legacy

FLAGS = $(CompilerFlags) -J$(ObjDir)

ObjsInCorrectLocation = $(addprefix $(ObjDir)/,$(Objs))
testObjsInCorrectLocation = $(addprefix $(ObjDir)/,$(testObjs))

# LIBS = -I$(MasterObjDir) -L$(MasterObjDir) -l$(MasterLibName)

# NetCDFlib = -I/usr/include -L/usr/lib/x86_64-linux-gnu -lnetcdff

all: $(ObjsInCorrectLocation)
	$(Compiler) $(FLAGS) $(StaticOpt) $^ -o $(Exec) # $(LIBS) $(NetCDFlib)

$(ObjDir)/%.o : $(CodeDir)/%.f
	$(Compiler) $(FLAGS) $(legacy) -c -o $@ $<

$(ObjDir)/%.o : $(CodeDir)/%.f90
	$(Compiler) $(FLAGS) -c -o $@ $<


test: $(testObjsInCorrectLocation)
	$(Compiler) $(FLAGS) $(StaticOpt) $^ -o testing # $(LIBS) $(NetCDFlib)

$(ObjDir)/%.o : $(CodeDir)/%.f
	$(Compiler) $(FLAGS) $(legacy) -c -o $@ $<

$(ObjDir)/%.o : $(CodeDir)/%.f90
	$(Compiler) $(FLAGS) -c -o $@ $< 


clean:
	rm -f $(ObjDir)/* $(Exec)





# cleanOutputs:
# 	rm -f buffer*.nc FinalSolVal*.csv ValFunCoefs*.nc fort.*

