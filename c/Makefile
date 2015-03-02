VERSION=1.5.3

#Compiler
CC = gcc -O3 -DVERSION='"$(VERSION)"' -g

#CC = gcc -O3 -DVERSION='"$(VERSION)"' -g

#compiler flags
# -g adds debug info to the executable file
# -Wall turns on most warnings from compiler
CFLAGS = -Wall

HTSLOC?=$(HTSLIB)

HTSTMP?=./htslib_tmp

#Define locations of header files
OPTINC?= -I$(HTSLOC)/
INCLUDES= -I./ $(OPTINC) -rdynamic

# define library paths in addition to /usr/lib
#   if I wanted to include libraries not in /usr/lib I'd specify
#   their path using -Lpath, something like:
LFLAGS?= -L$(HTSTMP)

# define any libraries to link into executable:
#   if I want to link in libraries (libx.so or libx.a) I use the -llibname
#   option, something like (this will link in libmylib.so and libm.so:
LIBS =-lhts -lpthread -lz -lm

# define the C source files
SRCS = ./bam_access.c ./bam_stats.c
#Define test sources
TEST_SRC=$(wildcard ./tests/*_tests.c)
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

#
# The following part of the makefile is generic; it can be used to
# build any executable just by changing the definitions above and by
# deleting dependencies appended to the file from 'make depend'
#

.PHONY: depend clean test make_htslib_tmp remove_htslib_tmp

.NOTPARALLEL: test

all: clean make_htslib_tmp $(BAM_STATS_TARGET) remove_htslib_tmp
	@echo  bam_stats compiled.

$(BAM_STATS_TARGET): $(OBJS)
	$(CC) $(CFLAGS) $(INCLUDES) -o $(BAM_STATS_TARGET) $(OBJS) $(LFLAGS) $(LIBS)

#Unit Tests
test: $(BAM_STATS_TARGET)
test: CFLAGS += $(INCLUDES) $(OBJS) $(LFLAGS) $(LIBS)
test: $(TESTS)
	sh ./tests/runtests.sh

#Unit tests with coverage
coverage: CFLAGS += --coverage
coverage: test

make_htslib_tmp:
	$(MD) $(HTSTMP)
	#Do some magic to ensure we compile BAM_STATS with the static libhts.a rather than libhts.so
	ln -s $(HTSLOC)/libhts.a $(HTSTMP)/libhts.a

remove_htslib_tmp:
	@echo remove tmp hts location
	-rm -rf $(HTSTMP)

copyscript:
	rsync -uE ./scripts/* ./bin/
	chmod a+x ../bin/bam_stats

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
	$(RM) ./*.o *~ $(BAM_STATS_TARGET) ./tests/tests_log $(TESTS) ./*.gcda ./*.gcov ./*.gcno *.gcda *.gcov *.gcno ./tests/*.gcda ./tests/*.gcov ./tests/*.gcno


depend: $(SRCS)
	makedepend $(INCLUDES) $^

# DO NOT DELETE THIS LINE -- make depend needs it
