#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright 2023 Ettus Research, a National Instruments Brand.
#
# SPDX-License-Identifier: GPL-3.0-or-later
#
#

from PyQt5 import Qt
from gnuradio import gr
import pmt


class StreamPushButton(gr.sync_block, Qt.QPushButton):
    """
    This block creates a variable push button that creates a message
    when clicked. The message will be formatted as a dictionary to pass
    to the RFNoC Replay block
    """

    def __init__(self, lbl, relBackColor, relFontColor,
                 mode, now, num_samps, message, chan=0, repeat=False):
        gr.sync_block.__init__(self, name="StreamPushButton",
                               in_sig=None, out_sig=None)
        Qt.QPushButton.__init__(self, lbl)

        self.lbl = lbl

        self.replayDict = {'stream_mode': mode}
        self.replayDict['stream_now'] = now
        self.replayDict['num_samps'] = num_samps
        self.replayDict['message'] = message
        self.replayDict['chan'] = chan

        self.sendDict = {'stream_cmd': self.replayDict}

        styleStr = ""
        if relBackColor != 'default':
            styleStr = "background-color: " + relBackColor + "; "

        if relFontColor:
            styleStr += "color: " + relFontColor + "; "

        self.setStyleSheet(styleStr)

        self.clicked[bool].connect(self.onBtnClicked)

        self.message_port_register_out(pmt.intern("pressed"))

    def set_mode(self, mode):
        self.replayDict['stream_cmd'] = mode

    def set_message(self, message):
        self.replayDict['message'] = message

    def set_chan(self, chan):
        self.replayDict['chan'] = chan

    def set_now(self, now):
        self.replayDict['stream_now'] = now

    def set_num_samps(self, num_samps):
        self.replayDict['num_samps'] = num_samps

    def onBtnClicked(self, pressed):
        self.message_port_pub(pmt.intern("pressed"), pmt.to_pmt(self.sendDict))
