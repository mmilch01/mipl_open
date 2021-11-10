include mipl.def
PROG	= dcm_sort/dcm_sort

include mlib3.mk
include dcmlib.def

CXXSRCS	= ${PROG}.cpp
OBJS	= ${CXXSRCS:.cpp=.o}
CXX = c++

ifeq ($(USE_4DFP),ON)
	LIBOBJS	= ${NMOBJS} ${TRXOBJS} ${DLIBOBJS} ${MLIB3OBJS} ${DCMLIBOBJS}
	CXXFLAGS = -O2 -fpermissive -Wno-deprecated  -I $(TRXDIR) -I $(MLIB3DIR) ${DCMTK_INCLUDES} -D_4DFP
else
	LIBOBJS = ${NMOBJS} ${DLIBOBJS} ${MLIB3OBJS} ${DCMLIBOBJS}
	CXXFLAGS = -O2 -fpermissive -Wno-deprecated  -I $(MLIB3DIR) ${DCMTK_INCLUDES}
endif

$(PROG): ${OBJS} ${LIBOBJS}
	$(CXX) ${DCMTK_CXXFLAGS} $(CXXFLAGS) ${DCMTK_LIBDIRS} -o $@ ${OBJS} ${LIBOBJS} ${DCMTK_LIBS}
.c.o:
	gcc -c $(CXXFLAGS) $< -o $@
.cpp.o:
	$(CXX) -c ${DCMTK_CXXFLAGS} $(CXXFLAGS) $< -o $@
clean:
	rm -f ${OBJS}
	rm -f ${PROG}
cleanall:
	rm -f ${OBJS}
	rm -f ${LIBOBJS}
	rm -f ${PROG}
release: ${PROG}
	chmod 771 ${PROG}
	cp -f ${PROG} ${RELEASE}
