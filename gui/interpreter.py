# In-process Python REPL for MFIX-GUI

from __future__ import print_function, absolute_import, unicode_literals, division

from qtpy import QtWidgets

from code import InteractiveConsole

# TODO:  colorize stderr, syntax highlighting, auto-complete, readline, clear window, search

class Interpreter(object):

    def init_interpreter(self):
        le = self.ui.lineedit_interpreter_input
        te = self.ui.textedit_interpreter_output
        le.editingFinished.connect(self.handle_interp_line)
        self.PS1 = '>>> '
        self.PS2 = '... '
        self.interpreter_prompt = self.PS1
        le.setText(self.interpreter_prompt)
        le.textEdited.connect(self.handle_interp_key)
        class Output:
            def __init__(self, textedit=te, err=False):
                self.textedit = te
                self.err = err
            def write(self, text):
                self.textedit.insertPlainText(text)
                #self.textedit.ensureCursorVisible()

        import sys
        sys.stdout = Output() # Will this grab all stdout/stderr?p
        sys.stderr = Output(err=True)
        self.interpreter = InteractiveConsole()
        banner = 'Python ' + sys.version + ' on ' + sys.platform + '\n'
        banner += 'MFIX-GUI version %s' % self.get_version() + '\n'
        banner += 'To access top-level MFIX-GUI object:\n'
        banner += '    from __main__ import gui\n'
        te.insertPlainText(banner)
        #self.textedit.ensureCursorVisible()

    def handle_interp_key(self):
        # Don't allow user to backspace over prompt
        #   TODO: readline/history support?
        le = self.ui.lineedit_interpreter_input
        text = le.text()
        if len(text) < 4:
            le.setText(self.interpreter_prompt)

    def handle_interp_line(self):
        le = self.ui.lineedit_interpreter_input
        te = self.ui.textedit_interpreter_output
        text = le.text().rstrip()
        if len(text) == 4 and self.interpreter_prompt == self.PS1: # No input
            return
        te.insertPlainText(text+"\n")
        text = text[4:] # Remove prompt
        result = self.interpreter.push(text)
        self.interpreter_prompt = self.PS2 if result else self.PS1
        le.setText(self.interpreter_prompt)
        te.ensureCursorVisible()
