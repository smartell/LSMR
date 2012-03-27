#Makefile for building and running the LSMR model on a Linux or Mac OSX
EXEC=LSMR
TPL=$(EXEC).tpl
DAT=simLSMR.dat

all: $(EXEC) $(EXEC).rep

$(EXEC): $(TPL)
	admb -s -g $(EXEC)

$(EXEC).rep: $(TPL) $(DAT)
	./$(EXEC) -ind $(DAT)

clean:
	rm -f $(EXEC) *.eva *.log *.rpt *.htp *.cor *.par *.r*
	rm -f *.rep *.bar *.psv *.std $(EXEC).cpp admodel.* variance