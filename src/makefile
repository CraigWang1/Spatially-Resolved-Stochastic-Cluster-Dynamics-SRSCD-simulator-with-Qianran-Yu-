# Specify which compiler to use:
CXX= mpicxx

# Use CXXFLAGS to specify additional compiler options to use:
CXXFLAGS = -Wall -std=c++17 -O3 -flto

# Define variables:
EXE_FILE= scdexe 
OBJ_FILES= main.o SCDWrapper.o Bundle.o CascadeDamage.o Damage.o Object.o OneLine.o rvgs.o cpdf.o BoundaryChange.o util.o
HEADER_FILES= SCDWrapper.h Bundle.h CascadeDamage.h Damage.h Object.h OneLine.h rvgs.h constants.h cpdf.h gnuplot-iostream.h gnuplot_i.h BoundaryChange.h util.h


# Build executable:
exe : $(EXE_FILE)
$(EXE_FILE) : $(OBJ_FILES) 
	$(CXX) -o $(EXE_FILE) $(OBJ_FILES) -lm -O3 -flto
	# chmod ugo+x $(EXE_FILE)	# turn on execute permissionss

# Rules:
.cpp.o:	$(HEADER_FILES)
	$(CXX) $(CXXFLAGS) -c  $*.cpp 

clean:
	rm $(EXE_FILE) $(OBJ_FILES)
