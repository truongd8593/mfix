.PHONY:  test test-gui tags

UIFILES := $(wildcard uifiles/*ui)
PY_UIFILES := $(patsubst %.ui,%.py,$(UIFILES))
PYUIC := $(shell { command -v pyuic4 || command -v pyuic5; } 2>/dev/null)

all: ${PY_UIFILES}

%.py: %.ui
	${PYUIC} -o $@ $<

clean:
	-rm -rf ${PY_UIFILES} dist __pycache__
	find . -name \*.pyc -delete

test2: all
	(Xvfb :99 &)
	DISPLAY=:99 MFIX_NO_VTK=1 python2 -m nose -v
	pkill Xvfb

test3: all
	(Xvfb :99 &)
	DISPLAY=:99 MFIX_NO_VTK=1 python3 -m nose -v
	pkill Xvfb

test: all
	(Xvfb :99 &)
	DISPLAY=:99 MFIX_NO_VTK=1 nosetests -v
	pkill Xvfb

test-gui:
	nosetests -v

tags:
	etags *.py {widgets,tests,tools}/*.py  # avoid uifiles

test-load-all:
	(Xvfb :99 &)
	find .. -name mfix.dat -print0 | DISPLAY=:99 MFIX_NO_VTK=1 xargs -0 -n 1 ./gui.py -lwarn -q &> /tmp/test-load-all.log & tail -f /tmp/test-load-all.log
	pkill Xvfb

pyinst: all
	pyinstaller gui.spec
