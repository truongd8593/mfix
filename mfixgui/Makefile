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
	MFIX_NO_VTK=1 xvfb-run python2 -m nose -v

test3: all
	MFIX_NO_VTK=1 xvfb-run python3 -m nose -v

test: all
	MFIX_NO_VTK=1 xvfb-run nosetests -v

test-gui:
	nosetests -v

tags:
	etags *.py {widgets,tests,tools}/*.py  # avoid uifiles

test-load-all:
	find .. -name mfix.dat -print0 | MFIX_NO_VTK=1 xvfb-run xargs -0 -n 1 ./run.sh -lwarning -q |& egrep --line-buf -v 'GL|Unsupported screen format' >& /tmp/test-load-all.log & tail -f /tmp/test-load-all.log

../mfix:
	cd ..; ./configure_mfix --python --smp
	make -C .. mfix

pyinst: all ../mfix
	pyinstaller gui.spec

sphinx-docs:
	cd doc/sphinx; sphinx-apidoc -o ../generated/source ../..
	cp doc/sphinx/index.rst doc/sphinx/conf.py doc/generated/source/
	cp doc/sphinx/Makefile doc/generated/
	make -C doc/generated html

docs: sphinx-docs
