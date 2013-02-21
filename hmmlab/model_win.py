from os.path import join, abspath, dirname
from gi.repository import GdkPixbuf
import gtklib

class ModelWindow(gtklib.ObjGetter):
    def __init__(self, model, main_win):
        path = join(dirname(abspath(__file__)), 'glade')
        gtklib.ObjGetter.__init__(self, join(path, "model_win.glade"), self.get_signals())
        self.model = model.model
        self.main_win = main_win
        self.window.show()

    def get_signals(self):
        signals = {"show" : self.show,
                "size_changed" : self.size_changed}
        return signals

    def show(self, widget):
        self.load()

    def load(self):
        alloc = self.image.get_allocation()
        path = self.model.create_image()
        pixbuf = GdkPixbuf.Pixbuf.new_from_file(path)
        width = pixbuf.get_width() if alloc.width > pixbuf.get_width() else alloc.width
        height = pixbuf.get_height() if alloc.height > pixbuf.get_height() else alloc.height
        pixbuf = GdkPixbuf.Pixbuf.new_from_file_at_scale(path,  width, height, False)
        self.image.set_from_pixbuf(pixbuf)

    def size_changed(self, widget, alloc):
        self.load()
