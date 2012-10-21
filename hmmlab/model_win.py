import gtklib
from os.path import join, abspath, dirname

class ModelWindow(gtklib.ObjGetter):
    def __init__(self, model, main_win):
        path = join(dirname(abspath(__file__)), 'glade')
        gtklib.ObjGetter.__init__(self, join(path, "model.glade"), self.get_signals())
        self.model = model
        self.main_win = main_win
        self.window.show()


