# Makefile

#############################################################
## modify to match your compiler defaults
CXX=g++

## modify to match your library paths
# VTK libraries
VTK := /cvmfs/soft.computecanada.ca/easybuild/software/2023/x86-64-v4/MPI/gcc12/openmpi4/paraview/5.13.1/

EXISTS := $(shell cmake --find-package -DNAME:STRING=VTK -DCOMPILER_ID:STRING=$(CXX) -DLANGUAGE:STRING=CXX -DMODE:STRING=EXIST &> /dev/null ; echo $$?)
ifeq ($(EXISTS),0)
VTK_CPPFLAGS := $(shell cmake --find-package -DNAME:STRING=VTK -DCOMPILER_ID:STRING=$(CXX) -DLANGUAGE:STRING=CXX -DMODE:STRING=COMPILE)
VTK_LDFLAGS := $(shell cmake --find-package -DNAME:STRING=VTK -DCOMPILER_ID:STRING=$(CXX) -DLANGUAGE:STRING=CXX -DMODE:STRING=LINK)
else
$(warning CMake not found or could not find VTK development headers. Using $$(VTK)=$(VTK) instead.)
VTK_CPPFLAGS := -I$(VTK)/include/paraview-5.13/ -I$(VTK)/include/paraview-5.13/vtk
#VTK_LDFLAGS := -L$(VTK)/lib64/vtk -L$(VTK)/lib//include/paraview-5.13/vtk -lvtkCommon -lvtkIO -lvtkFiltering -lvtkGraphics
VTK_LDFLAGS := -L$(VTK)/lib64/vtk -L$(VTK)/lib//include/paraview-5.13/vtk  \
  -lvtkCommonCore-pv5.13 \
  -lvtkCommonDataModel-pv5.13 \
  -lvtkFiltersCore-pv5.13 \
  -lvtkIOCore-pv5.13 \
  -lvtkIOXML-pv5.13 \
  -lvtkIOXMLParser-pv5.13 \
  -lvtkRenderingCore-pv5.13 \
  -lvtkRenderingOpenGL2-pv5.13 \
  -lvtkInteractionStyle-pv5.13 \
  -lvtksys-pv5.13 \
  -ldl -lpthread
endif

CPPFLAGS = -O3 $(VTK_CPPFLAGS) -Wall -Wno-deprecated
LIBS = $(VTK_LDFLAGS)
############################################################

TARGETS = mesh2vtu cell2vtu surf2vtu ugrid ugrid_pts

all: $(TARGETS)

%: %.o
	$(CXX) $(CPPFLAGS) -o $@ $< $(LIBS)

%.o: %.cxx
	$(CXX) $(CPPFLAGS) -c $< -o $@

clean:
	-rm -f $(TARGETS) $(addsuffix .o, $(TARGETS))
.PHONY: clean

