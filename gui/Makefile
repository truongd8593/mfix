UIFILES := $(wildcard uifiles/*ui)
PYFILES := $(patsubst %.ui,%.py,$(UIFILES))
PYUIC := $(shell pyuic5 --version >& /dev/null  && echo pyuic5 || pyuic4 --version >& /dev/null && echo pyuic4)

all: ${PYFILES}

%.py: %.ui
	${PYUIC} -o $@ $<

clean:
	rm -f ${PYFILES}

test:
	python -m unittest discover

nose:
	(Xvfb :99 &); \
	DISPLAY=:99 MFIX_NO_VTK=1 nosetests; \
	killall Xvfb
