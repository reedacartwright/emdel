AM_CPPFLAGS= -I$(top_srcdir)/include $(all_includes)

bin_PROGRAMS = emdel

emdel_CPPFLAGS = $(BOOST_CPPFLAGS) $(GSL_CFLAGS)
emdel_SOURCES = covar.cpp em.cpp emdel.cpp hist.cpp models.cpp seqdb.cpp \
	seqparser.cpp sample.cpp
emdel_LDFLAGS = $(all_libraries) $(BOOST_LDFLAGS) $(ALL_STATIC_LDFLAGS)
emdel_LDADD = ${LIBOJS} $(BOOST_PROGRAM_OPTIONS_LIB) $(BOOST_THREAD_LIB) $(GSL_LIBS)



#emdelwx_SOURCES = emdelwx.cpp
#emdelwx_LDFLAGS = $(all_libraries) 
#emdelwx_LDADD = ${LIBOJS} @LEXLIB@
#emdelwx_CPPFLAGS = ${CPPFLAGS} ${WX_CPPFLAGS@}
#emdelwx_CXXFLAGS = ${CXXFLAGS} ${WX_CXXFLAGS}
#emdelwx_CFLAGS = ${CFLAGS} ${WX_CFLAGS_ONLY}
#emdelwx_LIBS = ${LIBS} ${WX_LIBS}

noinst_HEADERS = ccvector.h covar.h covar_k2p.h em_k2p.h emdel.cmds emdel.h \
	hist.h hist_k2p.h seqdb.h table.h task.h models.h series.h invert_matrix.h \
	em.h big_prob.h sample.h sample_k2p.h seqparser.h xm.h
