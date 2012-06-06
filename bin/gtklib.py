from gi.repository import Gtk

class ObjGetter:
    def __init__(self, path, signals):
        self.builder = Gtk.Builder()
        self.builder.add_from_file(path)
        self.builder.connect_signals(signals)

    def __getattr__(self, name):
        return self.builder.get_object(name)
