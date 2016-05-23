.PHONY:  test test-gui

UIFILES := $(wildcard uifiles/*ui)
PY_UIFILES := $(patsubst %.ui,%.py,$(UIFILES))
PYUIC := $(shell { command -v pyuic4 || command -v pyuic5; } 2>/dev/null)

all: ${PY_UIFILES}

%.py: %.ui
	${PYUIC} -o $@ $<

clean:
	rm -f ${PY_UIFILES}
	find . -name \*.pyc -delete

test: all
	(Xvfb :99 &)
	DISPLAY=:99 MFIX_NO_VTK=1 nosetests -v
	pkill Xvfb

test-gui: all
	nosetests -v
