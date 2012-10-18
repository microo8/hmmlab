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
import math
from os.path import expanduser, join, exists, abspath, dirname
try:
    from hmmlablib import libhmm
    import gtklib
except ImportError:
    from hmmlab.hmmlablib import libhmm
    from hmmlab import gtklib

class DrawArea(gtklib.ObjGetter):
    '''Trieda drawingarea, ktora vykresluje jeden stream vo visualnom okne'''
    def __init__(self, stream_area):
        '''Vytvori visualne okno'''
        self.stream_area = stream_area
        path = join(dirname(abspath(__file__)), 'glade')
        gtklib.ObjGetter.__init__(self, join(path, 'drawarea.glade'), self.get_signals())

    def get_signals(self):
        signals = {'draw' : self.draw,
                   'set_wh' : self.set_wh}
        return signals

    def set_wh(self, widget, alloc):
        self.stream_area.set_wh(alloc.width, alloc.height)

    def draw(self, drawarea, cr):
        #clear
        cr.rectangle(0, 0, drawarea.get_allocation().width, drawarea.get_allocation().height)
        cr.set_source_rgb(0, 0, 0)
        cr.fill()
        cr.set_source_rgb (255, 200, 0)
        for pos in self.stream_area.pos_data:
            x = pos[0]
            y = pos[1]
            cr.move_to(x,y)
            cr.arc(x, y, 3, 0, 2*math.pi);
            cr.fill()
        cr.set_source_rgb (255, 0, 0)
        for pos in self.stream_area.pos_gaussians:
            x = pos[0]
            y = pos[1]
            cr.move_to(x,y)
            cr.arc(x, y, 3, 0, 2*math.pi);
            cr.fill()

