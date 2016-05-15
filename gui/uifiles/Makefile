UIFILES := $(wildcard *ui)
PYFILES := $(patsubst %.ui,%.py,$(UIFILES))
PYUIC := $(shell pyuic5 --version >& /dev/null  && echo pyuic5 || pyuic4 --version >& /dev/null && echo pyuic4)

all: ${PYFILES}

%.py: %.ui
	${PYUIC} -o $@ $<

clean:
	rm ${PYFILES}
