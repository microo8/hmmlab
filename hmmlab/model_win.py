'''
This file is part of HMMLab.

HMMLab is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

HMMLab is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with HMMLab.  If not, see <http://www.gnu.org/licenses/>.
'''

import os
from os.path import join, abspath, dirname
from gi.repository import Gtk, Gdk, GdkPixbuf
import gtklib

class ModelWindow(gtklib.ObjGetter):
    def __init__(self, model, main_win):
        path = join(dirname(abspath(__file__)), 'glade')
        gtklib.ObjGetter.__init__(self, join(path, "model_win.glade"), self.get_signals())
        self.model = model.model
        self.multiplier = 1
        self.image.override_background_color(Gtk.StateFlags.NORMAL, Gdk.RGBA(255,255,255))
        self.main_win = main_win
        self.window.show()

    def get_signals(self):
        signals = {"show" : self.show,
                "size_changed" : self.size_changed,
                "scroll" : self.scroll}
        return signals

    def scroll(self, widget, event):
        self.multiplier = max(0.1, self.multiplier - 0.1*event.delta_y)
        self.load()

    def show(self, widget):
        self.load()

    def load(self):
        alloc = self.image.get_allocation()
        path = self.model.create_image()
        pixbuf = GdkPixbuf.Pixbuf.new_from_file(path)
        width = pixbuf.get_width()
        height = pixbuf.get_height()
        if alloc.width < width:
            height *= alloc.width / width
            width = alloc.width
        if height > alloc.height:
            width = pixbuf.get_width()
            height = pixbuf.get_height()
            width *= alloc.height / height
            height = alloc.height
        pixbuf = GdkPixbuf.Pixbuf.new_from_file_at_scale(path, self.multiplier * width, self.multiplier * height, False)
        self.image.set_from_pixbuf(pixbuf)
        os.remove(path)

    def size_changed(self, widget, alloc):
        self.load()
