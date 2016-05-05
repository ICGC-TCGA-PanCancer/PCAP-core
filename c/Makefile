#VERSION=1.10.0
VERSION=$(shell perl -I../lib -MPCAP -e 'print PCAP->VERSION;')

#Compiler
CC = gcc -O3 -DVERSION='"$(VERSION)"' -g

#CC = gcc -O3 -DVERSION='"$(VERSION)"' -g

#compiler flags
# -g adds debug info to the executable file
# -Wall turns on most warnings from compiler
CFLAGS = -Wall

HTSLOC?=$(HTSLIB)
LIBBIGWIGLOC?=$(LIBBIGWIG)
HTSTMP?=./htslib_tmp
prefix=?/usr/local/

#Define locations of header files
OPTINC?= -I$(HTSLOC)/
INCLUDES= -I./ $(OPTINC) -rdynamic
CAT_OPTINC?= -I$(LIBBIGWIGLOC)
CAT_INCLUDES?= -I$(prefix)/include $(CAT_OPTINC)


# define library paths in addition to /usr/lib
#   if I wanted to include libraries not in /usr/lib I'd specify
#   their path using -Lpath, something like:
LFLAGS?= -L$(HTSTMP)
CAT_LFLAGS?= -L$(prefix)/lib -L$(LIBBIGWIGLOC)

# define any libraries to link into executable:
#   if I want to link in libraries (libx.so or libx.a) I use the -llibname
#   option, something like (this will link in libmylib.so and libm.so:
LIBS =-lhts -lpthread -lz -lm

# define the C source files
SRCS = ./bam_access.c ./bam_stats_output.c ./bam_stats_calcs.c
#Define test sources
TEST_SRC=$(wildcard ./c_tests/*_tests.c)
TESTS=$(patsubst %.c,%,$(TEST_SRC))

# define the C object files
#
# This uses Suffix Replacement within a macro:
#   $(name:string1=string2)
#         For each word in 'name' replace 'string1' with 'string2'
# Below we are replacing the suffix .c of all words in the macro SRCS
# with the .o suffix
#
OBJS = $(SRCS:.c=.o)

MD := mkdir

#Build target executable
BAM_STATS_TARGET=../bin/bam_stats
CAT_TARGET=../bin/bwcat
SQ_TARGET=../bin/reheadSQ

#
# The following part of the makefile is generic; it can be used to
# build any executable just by changing the definitions above and by
# deleting dependencies appended to the file from 'make depend'
#

.PHONY: depend clean test make_htslib_tmp remove_htslib_tmp

.NOTPARALLEL: test

all: clean make_htslib_tmp $(BAM_STATS_TARGET) test remove_htslib_tmp $(CAT_TARGET) $(SQ_TARGET)
	@echo  bam_stats, reheadSQ and bwcat compiled.

$(BAM_STATS_TARGET): $(OBJS)
	$(CC) $(CFLAGS) $(INCLUDES) -o $(BAM_STATS_TARGET) $(OBJS) $(LFLAGS) $(LIBS) ./bam_stats.c

$(CAT_TARGET):
	$(CC) $(CAT_INCLUDES) $(CFLAGS) ./bwcat.c $(CAT_LFLAGS) -lBigWig -lz -lcurl -lm -lgnutls -ltasn1 -lhogweed -lnettle -lgmp -lp11-kit -o $(CAT_TARGET)

$(SQ_TARGET):
	$(CC) $(CFLAGS) ./reheadSQ.c -o $(SQ_TARGET)


#Unit Tests
test: $(BAM_STATS_TARGET)
test: CFLAGS += $(INCLUDES) $(OBJS) $(LFLAGS) $(LIBS)
test: $(TESTS)
	sh ./c_tests/runtests.sh

#Unit tests with coverage
coverage: CFLAGS += --coverage
coverage: test

make_bigwigcat: $(CAT_TARGET)
	@echo $(CAT_TARGET) done

make_htslib_tmp:
	$(MD) $(HTSTMP)
	#Do some magic to ensure we compile BAM_STATS with the static libhts.a rather than libhts.so
	ln -s $(HTSLOC)/libhts.a $(HTSTMP)/libhts.a

remove_htslib_tmp:
	@echo remove tmp hts location
	-rm -rf $(HTSTMP)

copyscript:
	cp ./scripts/* ./bin/
	chmod a+x $(BAM_STATS_TARGET) $(CAT_TARGET) $(SQ_TARGET)

valgrind:
	VALGRIND="valgrind --log-file=/tmp/valgrind-%p.log" $(MAKE)


# this is a suffix replacement rule for building .o's from .c's
# it uses automatic variables $<: the name of the prerequisite of
# the rule(a .c file) and $@: the name of the target of the rule (a .o file)
# (see the gnu make manual section about automatic variables)
.c.o:
	$(CC) $(CFLAGS) $(INCLUDES) -c $<  -o $@

clean:
	@echo clean
	$(RM) ./*.o *~ $(BAM_STATS_TARGET) $(SQ_TARGET) $(CAT_TARGET) ./tests/tests_log $(TESTS) ./*.gcda ./*.gcov ./*.gcno *.gcda *.gcov *.gcno ./tests/*.gcda ./tests/*.gcov ./tests/*.gcno
	-rm -rf $(HTSTMP)


depend: $(SRCS)
	makedepend $(INCLUDES) $^

# DO NOT DELETE THIS LINE -- make depend needs it
