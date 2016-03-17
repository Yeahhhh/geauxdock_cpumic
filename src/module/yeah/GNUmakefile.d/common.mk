
MODULE := ./module

HEADPATH += -I$(MODULE)
HEADPATH += -I$(MODULE)/nvidia/cuda7
LIBPATH += -L$(MODULE)/lib
LINKFLAGS +=



