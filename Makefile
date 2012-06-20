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

# Set this variable if you want to use liblogger
LOGGER := 

# Set this variable if you want to use coverage
COVERAGE := 

ifdef LOGGER
 LOGGERINC = -I ../liblogger
 LOGGERLINK = -L../liblogger -llogger 
 LOGGERCFLAGS = -DUSE_LOGGER
endif

ifdef COVERAGE
	COVERAGECFLAGS = --coverage
	COVERAGELINK = -lgcov
endif

PLATFORM = $(shell uname)
CC = gcc

OBJ_DIR = obj
SRC_DIR = src
TEST_SRC_DIR = test_src
DEMO_SRC_DIR = demo_src
SOLVER_SRC_DIR = solver_src
BIN_DIR = bin
LIB_DIR = lib
DEPS_DIR = .deps
df = $(DEPS_DIR)/$(*F)

CFLAGS = -ggdb -g -pedantic -Wall -Wextra -Wswitch-default -Wswitch-enum -Wdeclaration-after-statement -Wmissing-declarations -Wmissing-include-dirs $(INCLUDES) $(LOGGERCFLAGS) $(COVERAGECFLAGS)
# this option not works for gcc 3.4.4
# -Wmissing-include-dirs
LIBCFLAGS = --std=c99
SOLVERCFLAGS = --std=gnu99

INCLUDES = -I . $(LOGGERINC)
LINKFLAGS = -L. -lspmatrix -lm $(LOGGERLINK) $(COVERAGELINK)
SOLVERLINKFLAGS = 

ifeq ($(PLATFORM),Linux)
LINKFLAGS += -lrt
endif

LIB_SOURCES = $(wildcard $(SRC_DIR)/*.c)
TEST_SOURCES = $(wildcard $(TEST_SRC_DIR)/*.c)	
DEMO_SOURCES = $(wildcard $(DEMO_SRC_DIR)/*.c)
SOLVER_SOURCES = $(wildcard $(SOLVER_SRC_DIR)/*.c)
LIB_OBJECTS = $(patsubst $(SRC_DIR)/%.c,$(OBJ_DIR)/%.o,$(LIB_SOURCES))
TEST_OBJECTS = $(patsubst $(TEST_SRC_DIR)/%.c,$(OBJ_DIR)/%.o,$(TEST_SOURCES))
DEMO_OBJECTS = $(patsubst $(DEMO_SRC_DIR)/%.c,$(OBJ_DIR)/%.o,$(DEMO_SOURCES))
SOLVER_OBJECTS = $(patsubst $(SOLVER_SRC_DIR)/%.c,$(OBJ_DIR)/%.o,$(SOLVER_SOURCES))

OBJECTS = $(LIB_OBJECTS) $(TEST_OBJECTS) $(DEMO_OBJECTS) $(SOLVER_OBJECTS)
DEPENDS = $(subst $(OBJ_DIR),$(DEPS_DIR),$(OBJECTS))

OUTPUT_TEST = $(BIN_DIR)/spmatrixtest
OUTPUT_LIB = $(LIB_DIR)/libspmatrix.a
FEM2D_DEMO = $(BIN_DIR)/demo_fem2d
OUTPUT_SOLVER = $(BIN_DIR)/solvertest
OUTPUT = $(OUTPUT_TEST) $(FEM2D_DEMO) $(OUTPUT_SOLVER)

# dependencies, based on article http://make.paulandlesley.org/autodep.html
# idea: generate dependencies, then replace "something.o : " with "something.o .deps/something.P : "
MAKEDEPEND = @gcc -MM $(CFLAGS) $(INCLUDES) -I src -o $*.d $<; sed 's/\($*\)\.o[ :]*/$(OBJ_DIR)\/\1.o $(DEPS_DIR)\/$*.P : /g' < $*.d > $*.P; rm $*.d; mv $*.P $(DEPS_DIR)

-include $(DEPENDS:%.o=%.P)

.DEFAULT_GOAL := all

# before starting the compilation be sure to create the objects directory
.PHONY:
all: $(OBJ_DIR) $(BIN_DIR) $(LIB_DIR) $(DEPS_DIR) $(OUTPUT)
	@echo "Build for $(PLATFORM) Done. See results in $(BIN_DIR) and $(LIB_DIR) directories"

# Rule for creation of the objects directory
$(OBJ_DIR):
	@mkdir -p $(OBJ_DIR)

$(BIN_DIR):
	@mkdir -p $(BIN_DIR)

$(LIB_DIR):
	@mkdir -p $(LIB_DIR)

$(DEPS_DIR):
	@mkdir -p $(DEPS_DIR)

# compile library sources
$(LIB_OBJECTS): $(OBJ_DIR)/%.o: $(SRC_DIR)/%.c
	$(MAKEDEPEND)
	$(CC) -c $(CFLAGS) $(LIBCFLAGS) $(DEFINES) $(INCLUDES) -c $< -o $@

# compile test sources
$(TEST_OBJECTS): $(OBJ_DIR)/%.o:$(TEST_SRC_DIR)/%.c
	$(MAKEDEPEND)
	$(CC) $(CFLAGS) $(LIBCFLAGS) $(DEFINES) $(INCLUDES) -I $(SRC_DIR) -c $< -o $@

# compile demo sources
$(DEMO_OBJECTS): $(OBJ_DIR)/%.o:$(DEMO_SRC_DIR)/%.c
	$(MAKEDEPEND)
	$(CC) $(CFLAGS) $(LIBCFLAGS) $(DEFINES) $(INCLUDES) -I $(SRC_DIR) -c $< -o $@

# compile solver test sources
$(SOLVER_OBJECTS): $(OBJ_DIR)/%.o:$(SOLVER_SRC_DIR)/%.c
	$(MAKEDEPEND)
	$(CC) $(CFLAGS) $(SOLVERCFLAGS) $(DEFINES) $(INCLUDES) -I $(SRC_DIR) -c $< -o $@


# link binaries
$(OUTPUT_TEST): $(OUTPUT_LIB) $(BIN_DIR) $(TEST_OBJECTS)
	$(CC) $(TEST_OBJECTS) -o $(OUTPUT_TEST) $(LINKFLAGS) -L $(LIB_DIR)

$(OUTPUT_SOLVER): $(OUTPUT_LIB) $(BIN_DIR) $(SOLVER_OBJECTS)
	$(CC) $(SOLVER_OBJECTS) -o $(OUTPUT_SOLVER) $(LINKFLAGS) $(SOLVERLINKFLAGS) -L $(LIB_DIR) 

$(FEM2D_DEMO): $(OUTPUT_LIB) $(BIN_DIR) $(DEMO_OBJECTS)
	$(CC) $(DEMO_OBJECTS) -o $(FEM2D_DEMO) $(LINKFLAGS) -L $(LIB_DIR)

# link the library
$(OUTPUT_LIB): $(LIB_OBJECTS) $(LIB_DIR)
	$(RM) -f $(OUTPUT_LIB)
	$(AR) cr $(OUTPUT_LIB) $(LIB_OBJECTS)
	ranlib $(OUTPUT_LIB)

lint:
	splint $(LIB_SOURCES)

test: $(OUTPUT_TEST)
	./$(OUTPUT_TEST)

.PHONY : clean
clean :
	@rm -f $(LIB_OBJECTS) $(TEST_OBJECTS) $(DEMO_OBJECTS) $(SOLVER_OBJECTS)
	@rm -f $(OUTPUT_TEST) $(FEM2D_DEMO) $(OUTPUT_LIB) $(OUTPUT_SOLVER)
	@rmdir $(BIN_DIR) $(LIB_DIR) $(OBJ_DIR)
	@rm -fr $(DEPS_DIR)

check-syntax: 
	gcc -o nul -S ${CHK_SOURCES} 

