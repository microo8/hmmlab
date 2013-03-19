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
    def __init__(self, main_win, modelset):
        '''Vytvori visualne okno'''
        path = join(dirname(abspath(__file__)), 'glade')
        gtklib.ObjGetter.__init__(self, join(path, 'visual_win.glade'), self.get_signals())
        self.config = configparser.ConfigParser()
        self.config.read(expanduser('~/.config/hmmlab.conf'))
        self.window.set_default_size(int(self.config['visualwindow']['width']), int(self.config['visualwindow']['height']))
        self.main_win = main_win
        self.streams = []
        self.toggled = True
        self.set_modelset(modelset)
        self.toggle(self.togglebutton1, True)
        if self.modelset is not None:
            self.adjustment1.set_upper(self.modelset.streams_distribution[0])
            self.adjustment2.set_upper(self.modelset.streams_distribution[0])

    def set_modelset(self, modelset):
        assert(self.modelset is None)
        self.modelset = modelset
        if self.modelset is not None:
            self.adjustment1.set_upper(self.modelset.streams_distribution[0])
            self.adjustment2.set_upper(self.modelset.streams_distribution[0])
            for strarea in self.modelset.stream_areas:
                da = DrawArea(self.main_win, strarea)
                da.adjustment1 = self.adjustment1
                da.adjustment2 = self.adjustment2
                self.vbox.pack_start(da.eventbox, True, True, 3)
                self.streams.append(da)
            self.window.show()
            self.refresh()

    def get_signals(self):
        signals = {'destroy': self.destroy,
                   'graphviz_toggled' : self.graphviz_toggled,
                   'pca_toggled' : self.pca_toggled,
                   'prob_toggled' : self.prob_toggled,
                   'prierez_toggled' : self.prierez_toggled,
                   'dim_changed' : self.dim_changed}
        return signals

    def destroy(self, *args):
        self.window.destroy()

    def __del__(self):
        self.destroy()

    def refresh(self, wh=False):
        for stream in self.streams:
            stream.drawarea.queue_draw()

    def toggle(self, button, state):
        self.toggled = False
        button.set_active(state)
        self.toggled = True

    def dim_changed(self, adj):
        for stream in self.streams:
            stream.drawarea.queue_draw()

    def graphviz_toggled(self, button=None):
        if self.toggled:
            if self.togglebutton1.get_active():
                self.toggle(self.togglebutton2, False)
                self.toggle(self.togglebutton3, False)
                self.toggle(self.togglebutton4, False)
                self.witch_dims.hide()
                for stream in self.streams:
                    stream.state = 'graphviz'
                    stream.drawarea.queue_draw()

    def pca_toggled(self, button=None):
        if self.toggled:
            if self.togglebutton2.get_active():
                self.toggle(self.togglebutton1, False)
                self.toggle(self.togglebutton3, False)
                self.toggle(self.togglebutton4, False)
                self.witch_dims.hide()
                for stream in self.streams:
                    stream.state = 'pca'
                    stream.drawarea.queue_draw()

    def prob_toggled(self, button=None):
        if self.toggled:
            if self.togglebutton3.get_active():
                self.toggle(self.togglebutton1, False)
                self.toggle(self.togglebutton2, False)
                self.toggle(self.togglebutton4, False)
                self.witch_dims.hide()
                for stream in self.streams:
                    stream.state = 'prob'
                    stream.drawarea.queue_draw()

    def prierez_toggled(self, button=None):
        if self.toggled:
            if self.togglebutton4.get_active():
                self.toggle(self.togglebutton1, False)
                self.toggle(self.togglebutton2, False)
                self.toggle(self.togglebutton3, False)
                self.witch_dims.show()
                for stream in self.streams:
                    stream.state = 'prierez'
                    stream.drawarea.queue_draw()


