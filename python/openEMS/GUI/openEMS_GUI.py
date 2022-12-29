# -*- coding: utf-8 -*-
"""
Created on Wed Jan  1 16:46:52 2020

@author: thorsten.liebig@gmx.de
"""

PYSIDE=False

if PYSIDE:
    from PySide2.QtGui import *
    from PySide2.QtCore import *
    from PySide2.QtWidgets import *
    from PySide2.QtCore import Signal
    pyqtSignal = Signal
else:
    from PyQt5.QtGui import *
    from PyQt5.QtCore import *
    from PyQt5.QtWidgets import *

from threading import Thread
import os, sys

WARN_COLOR = 'orange'
ERR_COLOR = 'red'

def createApp(argv=None):
    if argv is None:
        argv = sys.argv

    app = QApplication(argv)
    return app

class logger:
    def __init__(self, log_cb):
        self.log_cb = log_cb

    def write(self, msg):
        self.log_cb(msg)

    def flush(self):
        pass

import threading
import time

# modified and fixed from:
# https://stackoverflow.com/questions/24277488/in-python-how-to-capture-the-stdout-from-a-c-shared-library-to-a-variable
class OutputGrabber():
    """
    Class used to grab standard output or another stream.
    """
    escape_char = "\b"

    def __init__(self, stream, log_cb):
        self.origstream = stream
        if self.origstream is None:
            self.origstream = sys.stdout
        self.origstreamfd = self.origstream.fileno()
        self.capturedtext = ""
        # Create a pipe so the stream can be captured:
        self.pipe_out, self.pipe_in = os.pipe()
        self.log_cb = log_cb

    def __enter__(self):
        self.start()
        return self

    def __exit__(self, type, value, traceback):
        self.stop()

    def start(self):
        """
        Start capturing the stream data.
        """
        self.capturedtext = ""
        # Save a copy of the stream:
        self.streamfd = os.dup(self.origstreamfd)
        # Replace the original stream with our write pipe:
        os.dup2(self.pipe_in, self.origstreamfd)
        # Start thread that will read the stream:
        self.workerThread = threading.Thread(target=self.readOutput)
        self.workerThread.start()
        # Make sure that the thread is running and os.read() has executed:
        time.sleep(0.01)

    def stop(self):
        """
        Stop capturing the stream data and save the text in `capturedtext`.
        """
        # give some time for the last output to be processed
        time.sleep(0.1)
        # Print the escape character to make the readOutput method stop:
        os.write(self.pipe_in,self.escape_char.encode(self.origstream.encoding))
        # wait until the thread finishes so we are sure that
        # we have until the last character:
        self.workerThread.join()
        # Close the pipe:
        os.close(self.pipe_in)
        os.close(self.pipe_out)
        # Restore the original stream:
        os.dup2(self.streamfd, self.origstreamfd)
        # Close the duplicate stream:
        os.close(self.streamfd)

    def readOutput(self):
        """
        Read the stream data (one byte at a time)
        and save the text in `capturedtext`.
        Sent it to the log callback in case of a newline
        """
        while True:
            char = os.read(self.pipe_out, 1).decode(self.origstream.encoding)
            if not char or self.escape_char in char:
                if len(self.capturedtext)>0:
                    self.log_cb(self.capturedtext)
                break
            self.capturedtext += char
            if char == '\n':
                self.log_cb(self.capturedtext)
                self.capturedtext = ''

class openEMS_Thread(QThread):
    def __init__(self, FDTD, sim_path, log_cb, run_kw={}):
        super(openEMS_Thread, self).__init__()
        self.FDTD     = FDTD
        self.sim_path = sim_path
        self.run_kw   = run_kw
        self.log_cb   = log_cb
        self.std_out_grab = None
        self.std_err_grab = None

    def run(self):
        self.std_out_grab = OutputGrabber(sys.stdout, self.log_cb)
        self.std_out_grab.start()

        self.std_err_grab = OutputGrabber(sys.stderr, lambda msg: self.log_cb(msg, ERR_COLOR))
        self.std_err_grab.start()

        self.FDTD.Run(self.sim_path, **self.run_kw)

        self.std_out_grab.stop()
        self.std_err_grab.stop()


class openEMS_CTRL(QDialog):
    sig_log_msg  = pyqtSignal(str, str)
    def __init__(self, parent, FDTD, sim_path, auto_close=True, **run_kw):
        super(openEMS_CTRL, self).__init__(parent)
        self.setWindowTitle('openEMS RUN Control')
        self.FDTD     = FDTD
        self.auto_close = auto_close
        self.was_aborted = False

        self.main_vlay = QVBoxLayout()
        self.log_win   = QTextEdit()
        self.log_win.setFont(QFont('DejaVu Sans Mono'))
        self.log_win.setReadOnly(True)
        self.log_win.setLineWrapMode(QTextEdit.NoWrap)
        self.sig_log_msg.connect(self._addLogLine)

        self.main_vlay.addWidget(self.log_win, stretch=1)

        self.buttonBox = QDialogButtonBox()
        self.buttonBox.accepted.connect(self.accept)
        self.buttonBox.rejected.connect(self.reject)
#        self.ok_button     = self.buttonBox.addButton(QDialogButtonBox.Ok)
        self.cancel_button = self.buttonBox.addButton(QDialogButtonBox.Cancel)
        self.cancel_button.setText('Abort')
        self.main_vlay.addWidget(self.buttonBox, stretch=0)
        self.setLayout(self.main_vlay)

        self.fdtd_thread = openEMS_Thread(FDTD, sim_path, self.addLogLine, run_kw)
        self.fdtd_thread.finished.connect(self.simFinished)
        self.fdtd_thread.start()

        self.resize(900, 600)

    def simFinished(self):
        if self.auto_close:
            self.close()
            return
        self.cancel_button.setEnabled(True)
        self.cancel_button.setText('Close')

    def _addLogLine(self, msg, color):
        # internal not thread save add log line
        textformat = QTextCharFormat()
        textformat.setForeground(QBrush(QColor(color)))
        # remember if scrollbar was at the bottom
        bottom = self.log_win.verticalScrollBar().sliderPosition() > (self.log_win.verticalScrollBar().maximum()-20)
        tc = QTextCursor(self.log_win.textCursor())
        tc.movePosition(QTextCursor.End)
        tc.insertText(msg, textformat)
        # if scrollbar was at the bottom, keep it there
        if bottom:
            self.log_win.verticalScrollBar().setValue(self.log_win.verticalScrollBar().maximum())

    def addLogLine(self, msg, color=None):
        if color is None:
            l_msg = msg.lower()
            if 'error' in l_msg:
                color = ERR_COLOR
            elif 'warning' in l_msg:
                color = WARN_COLOR
        self.sig_log_msg.emit(msg, color)

    def wasAborted(self):
        return self.was_aborted

    def accept(self):
        if not self.fdtd_thread.isFinished():
            return
        super(openEMS_CTRL, self).accept()

    def reject(self):
        if not self.fdtd_thread.isFinished():
            ret =  QMessageBox.question(self, 'Abort Simulation?', 'Abort the current simulation?', QMessageBox.Yes | QMessageBox.No)
            if ret != QMessageBox.Yes:
                return
            self.was_aborted = True
            self.FDTD.SetAbort(True)
            self.cancel_button.setEnabled(False)
            return
        super(openEMS_CTRL, self).reject()

def runOpenEMS(FDTD, sim_path, use_GUI=True, auto_close=True, **kw):
    if not use_GUI:
        return FDTD.Run(sim_path, **kw)

    app = createApp()

    gui = openEMS_CTRL(None, FDTD, sim_path, auto_close=auto_close, **kw)

    gui.exec()

    return not gui.wasAborted()
