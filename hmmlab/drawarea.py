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
try:
    from hmmlablib import libhmm
    import gtklib
except ImportError:
    from hmmlab.hmmlablib import libhmm
    from hmmlab import gtklib

class DrawArea(gtklib.ObjGetter):
    '''Trieda drawingarea, ktora vykresluje jeden stream vo visualnom okne'''
    def __init__(self, i, modelset=None):
        '''Vytvori visualne okno'''
        path = join(dirname(abspath(__file__)), 'glade')
        gtklib.ObjGetter.__init__(self, join(path, 'drawarea.glade'), self.get_signals())
        self.i = i
        self.modelset = modelset

    def get_signals(self):
        signals = {'draw' : self.draw}
        return signals

    def draw(self, drawarea, cr):
        #clear
        cr.rectangle(0, 0, drawarea.get_allocation().width, drawarea.get_allocation().height)
        cr.set_source_rgb(255, 255, 255)
        cr.fill()
