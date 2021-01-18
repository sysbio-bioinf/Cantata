SRCDIR = src
OBJDIR = obj
BINDIR = bin
OBJ = $(OBJDIR)/booleanformula.o $(OBJDIR)/constraintviolation.o $(OBJDIR)/formulaparser.o $(OBJDIR)/gp.o $(OBJDIR)/helpers.o $(OBJDIR)/individual.o $(OBJDIR)/main.o $(OBJDIR)/mersennetwister.o $(OBJDIR)/networkconstraint.o 

CC = g++
FLAGS = -Wall -O3 -g -c
LNFLAGS = 
OUT = cantata
INCLUDE = /usr/local/boost

all: $(OBJ)
	@echo 'Linking object files to $(BINDIR)/$(OUT)...'
	@if [ ! -d $(BINDIR) ]; then mkdir -p $(BINDIR); fi
	@$(CC) $(LNFLAGS) -o "$(BINDIR)/$(OUT)" $^
	@echo 'Build finished!'

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	@echo 'Building file: $@...'
	@if [ ! -d $(OBJDIR) ]; then mkdir -p $(OBJDIR); fi
	@$(CC) -I$(INCLUDE) $(FLAGS) $< -o $@

.PHONY: clean all
clean:
	@echo 'Removing binaries...'
	@rm -f $(OBJ)
	@if [ -d $(OBJDIR) ]; then rmdir -p --ignore-fail-on-non-empty $(OBJDIR); fi
	@rm -f $(BINDIR)/$(OUT)*
	@if [ -d $(BINDIR) ]; then rmdir -p $(BINDIR); fi
	@echo 'Cleaning finished!'
