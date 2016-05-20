.PHONY:  test nose

UIFILES := $(wildcard uifiles/*ui)
PY_UIFILES := $(patsubst %.ui,%.py,$(UIFILES))
PYUIC := $(shell pyuic5 --version >& /dev/null  && echo pyuic5 || pyuic4 --version >& /dev/null && echo pyuic4)

all: ${PY_UIFILES}

%.py: %.ui
	${PYUIC} -o $@ $<

clean:
	rm -f ${PY_UIFILES}
	find . -name \*.pyc -delete

test:   all
	python -m unittest discover

nose:   all
	(Xvfb :99 &); \
	DISPLAY=:99 MFIX_NO_VTK=1 nosetests -v; \
	killall Xvfb
