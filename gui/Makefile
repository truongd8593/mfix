.PHONY:  test test-gui tags

UIFILES := $(wildcard uifiles/*ui)
PY_UIFILES := $(patsubst %.ui,%.py,$(UIFILES))
PYUIC := $(shell { command -v pyuic5 || command -v pyuic4; } 2>/dev/null)

all: ${PY_UIFILES}

%.py: %.ui
	${PYUIC} -o $@ $<

clean:
	-rm -rf ${PY_UIFILES} dist __pycache__
	rm -rf doc/generated
	find . -name \*.pyc -delete
	make -C .. clean

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
	find .. -name mfix.dat -print0 | DISPLAY=:99 MFIX_NO_VTK=1 xargs -0 -n 1 ./gui.py -lwarn -q |& egrep --line-buf -v 'GL|Unsupported screen format' >& /tmp/test-load-all.log & tail -f /tmp/test-load-all.log
	pkill Xvfb

../mfix:
	cd ..; ./configure_mfix --python --smp
	make -C .. mfix

pyinst: all ../mfix
	pyinstaller gui.spec

sphinx-docs:
	cd doc/sphinx; sphinx-apidoc -o ../generated/source ../..
	cp doc/sphinx/index.rst doc/sphinx/conf.py doc/generated/source
	make -C doc/sphinx html

docs: sphinx-docs
