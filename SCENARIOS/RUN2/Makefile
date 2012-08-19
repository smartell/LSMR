## Makefile for running iscam
## Targets: 
##		all:   -copy executable and run the model with DAT & ARG
##		run:   -copy executable and force a run
##		mcmc:  -copy executable and run in mcmc mode and mceval
##		retro: -copy executable and run  retrospective analysis
EXEC=lsmr
prefix=../../ADMB/srcLSMR
DAT=RUN2.dat
CTL=RUN2
ARG=
MCFLAG=-mcmc 10000 -mcsave 100 -nosdmcmc
NR=4

ifdef DEBUG
  DIST=$(prefix)/debug/$(EXEC)
else
  DIST=$(prefix)/release/$(EXEC)
endif


all: $(EXEC) $(EXEC).rep dust

$(EXEC): $(DIST)
	cp $(DIST) $@

$(EXEC).rep: $(DIST) $(CTL).ctl 
	./$(EXEC) -ind $(DAT) $(ARG)


run:  $(EXEC)
	./$(EXEC) -ind $(DAT) $(ARG)
	make dust

mcmc: $(EXEC) $(EXEC).psv
	./$(EXEC) -ind $(DAT) -mceval

$(EXEC).psv: $(CTL).ctl
	./$(EXEC) -ind $(DAT) $(MCFLAG) $(ARG)

mceval: $(EXEC)
	cp $(CTL).psv $(EXEC).psv
	./$(EXEC) -ind $(DAT) -mceval

retro: $(EXEC) $(EXEC).ret1

$(EXEC).ret1:
	@echo $(RUNRETRO) | R --vanilla --slave

RUNRETRO = 'args = paste("-retro",c(1:$(NR),0)); \
            sapply(args,\
            function(a){ cmd=paste("./$(EXEC)","-ind $(DAT)",a);\
                        system(cmd)})'

dust:
	-rm -f admodel.* $(EXEC).[bel]* $(EXEC).*[0123456789]*

clean: 
	-rm -f $(EXEC).* admodel.* variance eigv.rpt fmin.log $(EXEC) variance