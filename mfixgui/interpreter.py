# In-process Python REPL for MFIX-GUI

from __future__ import print_function, absolute_import, unicode_literals, division

from code import InteractiveConsole

# TODO:  colorize stderr, syntax highlighting, auto-complete, readline, clear window, search

class Interpreter(object):

    def init_interpreter(self):
        import sys
        self.interp_setup_done = False
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
                text = text.rstrip()
                if text:
                    self.textedit.appendPlainText(text)

        self.stdout = Output()
        self.stderr = Output(err=True)
        self.interpreter = InteractiveConsole()
        banner = 'Python ' + sys.version + ' on ' + sys.platform + '\n'
        banner += 'MFIX-GUI version %s' % self.get_version() + '\n'
        te.insertPlainText(banner)

    def setup_interpreter(self):
        if self.interp_setup_done:
            return

        self.interp_setup_done = True
        for line in ("from __main__ import gui as g", "p = g.project"):
            self.ui.lineedit_interpreter_input.setText('>>> ' + line)
            self.handle_interp_line()

    def capture_output(self, enable):
        import sys
        if enable:
            sys.stdout = self.stdout
            sys.stderr = self.stderr
        else:
            sys.stdout = sys.__stdout__
            sys.stderr = sys.__stderr__

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
        if len(text) == 3 and self.interpreter_prompt == self.PS1: # No input
            return
        te.appendPlainText(text)
        text = text[4:] # Remove prompt
        result = self.interpreter.push(text)
        self.interpreter_prompt = self.PS2 if result else self.PS1
        le.setText(self.interpreter_prompt)
