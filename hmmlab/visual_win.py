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

from os.path import expanduser, join, exists, abspath, dirname
import configparser
import cairo
try:
    from hmmlablib import libhmm
    from drawarea import DrawArea
    import gtklib
except ImportError:
    from hmmlab.hmmlablib import libhmm
    from hmmlab.drawarea import DrawArea
    from hmmlab import gtklib


class VisualWindow(gtklib.ObjGetter):
    '''Trieda visualneho okna'''
    def __init__(self, modelset):
        '''Vytvori visualne okno'''
        path = join(dirname(abspath(__file__)), 'glade')
        gtklib.ObjGetter.__init__(self, join(path, 'visual_win.glade'), self.get_signals())
        self.config = configparser.ConfigParser()
        self.config.read(expanduser('~/.config/hmmlab.conf'))
        self.window.set_default_size(int(self.config['visualwindow']['width']), int(self.config['visualwindow']['height']))
        self.streams = []
        self.set_modelset(modelset)

    def set_modelset(self, modelset):
        assert(self.modelset is None)
        self.modelset = modelset
        if self.modelset is not None:
            for i in range(self.modelset.streams_size):
                da = DrawArea(i, self.modelset)
                self.vbox.pack_start(da.eventbox, True, True, 3)
                self.streams.append(da)
            self.window.show()
            self.refresh()

    def get_signals(self):
        signals = {'destroy': self.destroy}
        return signals

    def destroy(self, *args):
        self.window.destroy()

    def __del__(self):
        print('asd')
        self.destroy()

    def refresh(self):
        for stream in self.streams:
            stream.drawarea.queue_draw()
