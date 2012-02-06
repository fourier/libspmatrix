# -*- Mode: makefile; -*-

# Copyright (C) 2011,2012 Alexey Veretennikov (alexey dot veretennikov at gmail.com)
# 
#	This file is part of libspmatrix.
#
# libspmatrix is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# libspmatrix is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with libspmatrix.  If not, see <http://www.gnu.org/licenses/>.


CC = gcc

CFLAGS = -ggdb -g --std=c99 -pedantic -Wall -Wextra -Wswitch-default -Wswitch-enum -Wdeclaration-after-statement -Wmissing-declarations $(INCLUDES)
# this option not works for gcc 3.4.4
#-Wmissing-include-dirs
INCLUDES = -I . -I ../liblogger
LINKFLAGS = -L. -lspmatrix -lm -L../liblogger -llogger

OUTPUT_SRC = main.c
FEM2D_SRC = main_fem2d.c
SOURCES := $(wildcard *.c)
HEADERS := $(wildcard *.h)
OBJECTS := $(patsubst %.c,%.o,$(SOURCES))
OBJECTS_LIB := $(filter-out $(patsubst %.c,%.o,$(OUTPUT_SRC)),$(OBJECTS))
OUTPUT_TEST = spmatrixtest
OUTPUT_LIB = libspmatrix.a
FEM2D_DEMO = demo_fem2d
OUTPUT = $(OUTPUT_TEST) $(FEM2D_DEMO)

all: $(OUTPUT)
	@echo "Done"

%.o : %.c %.h
	$(CC) -c $(CFLAGS) $(DEFINES) $(INCLUDES) $< -o $@


$(OUTPUT_TEST): $(OUTPUT_LIB) 
	$(CC) $(patsubst %.c,%.o,$(OUTPUT_SRC)) -o $(OUTPUT_TEST) $(LINKFLAGS)


$(FEM2D_DEMO): $(OUTPUT_LIB) 
	$(CC) $(patsubst %.c,%.o,$(FEM2D_SRC)) -o $(FEM2D_DEMO) $(LINKFLAGS)


$(OUTPUT_LIB): $(OBJECTS)
	$(RM) -f $(OUTPUT_LIB)
	$(AR) cr $(OUTPUT_LIB) $(OBJECTS_LIB)
	ranlib $(OUTPUT_LIB)

lint:
	splint *.c


.PHONY : clean
clean :
	rm $(OBJECTS) $(OUTPUT_TEST) $(FEM2D_DEMO) $(OUTPUT_LIB)

check-syntax: 
	gcc -o nul -S ${CHK_SOURCES} 

