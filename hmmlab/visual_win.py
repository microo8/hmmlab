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

from os.path import expanduser, join, exists
import cairo
try:
    from hmmlablib import libhmm
    import gtklib
except ImportError:
    from hmmlab.hmmlablib import libhmm
    from hmmlab import gtklib


class VisualWindow(gtklib.ObjGetter):
    '''Trieda visualneho okna'''
    def __init__(self, modelset=None):
        '''Vytvori visualne okno'''
        path = join(os.path.dirname(os.path.abspath(__file__)), 'glade')
        gtklib.ObjGetter.__init__(self, join(path, 'visual_win.glade'), self.get_signals())]
        self.modelset = modelset
        self.window.show()

    def get_signals(self):
        signals = {}
        return signals

    def destroy(self):
        self.window.destroy()

    def __del__(self):
        self.destroy()

    def draw(self, drawarea, cr):
        if self.modelset is not None:
            #TODO
            pass
