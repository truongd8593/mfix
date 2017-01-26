.PHONY:  test test-gui tags

UIFILES := $(wildcard uifiles/*ui)
PY_UIFILES := $(patsubst %.ui,%.py,$(UIFILES))
PYUIC := $(shell { command -v pyuic5 || command -v pyuic4; } 2>/dev/null)

all: ${PY_UIFILES}

%.py: %.ui
	${PYUIC} -o $@ $<

clean:
	rm -rf ${PY_UIFILES}
	rm -rf dist
	rm -rf doc/generated
	rm -rf __pycache__ */__pycache__
	find . -name \*.pyc -delete
	find .. -name \*.stl.cgs -delete
	find .. -name particle_input.dat.cgs -delete
	find .. -name geometry.stl.original -delete

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

test-gui:   # test with GUI enabled
	nosetests -v

tags:
	etags *.py {widgets,tests,tools}/*.py  # avoid uifiles

test-load-all: clean test-load-all-2 test-load-all-3

test-load-all-2: clean
	find .. -name mfix.dat -print0 | MFIX_NO_VTK=1 xargs -0 -n 1 ./run2.sh -linfo -t
test-load-all-3: clean
	find .. -name mfix.dat -print0 | MFIX_NO_VTK=1 xargs -0 -n 1 ./run3.sh -linfo -t

../mfix:
	cd ..; ./configure_mfix --python --smp
	make -C .. mfix mfixsolver.so

pyinst: all ../mfix
	pyinstaller gui.spec

sphinx-docs:
	cd doc/sphinx; sphinx-apidoc -o ../generated/source ../..
	cp doc/sphinx/index.rst doc/sphinx/conf.py doc/generated/source/
	cp doc/sphinx/Makefile doc/generated/
	make -C doc/generated html

docs: sphinx-docs
