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
        self.selected_gaussian_index = None
        self.state = 'graphviz'
        self.draw_functions = {'graphviz' : '',
                               'pca' : '_pca'}

    def get_signals(self):
        signals = {'draw' : self.draw,
                   'set_wh' : self.set_wh,
                   'press' : self.press}
        return signals

    def set_wh(self, widget, alloc):
        self.stream_area.set_wh(alloc.width, alloc.height)

    def draw(self, drawarea, cr):
        #clear
        cr.rectangle(0, 0, drawarea.get_allocation().width, drawarea.get_allocation().height)
        cr.set_source_rgb(0, 0, 0)
        cr.fill()
        pos_list = getattr(self.stream_area, 'pos_data' + self.draw_functions[self.state])
        for i, pos in enumerate(pos_list):
            x = pos[0]
            y = pos[1]
            if self.selected_gaussian_index is not None and i in self.stream_area.selected_gaussians[self.selected_gaussian_index].my_data:
                    cr.set_source_rgb(0, 100, 200)
            else:
                cr.set_source_rgb (255, 200, 0)
            cr.move_to(x,y)
            cr.arc(x, y, 3, 0, 2*math.pi);
            cr.fill()
        pos_list = getattr(self.stream_area, 'pos_gaussians' + self.draw_functions[self.state])
        for i, pos in enumerate(pos_list):
            x = pos[0]
            y = pos[1]
            if self.selected_gaussian_index == i:
                cr.set_source_rgb(0, 0, 255)
            else:
                cr.set_source_rgb (255, 0, 0)
            cr.move_to(x,y)
            cr.arc(x, y, 3, 0, 2*math.pi);
            cr.fill()
            if self.state == 'pca':
                ew = self.stream_area.pos_gaussians_var_pca[i][0]
                eh = self.stream_area.pos_gaussians_var_pca[i][1]
                gtklib.cairo_ellipse(cr, x - ew / 2., y - eh / 2., ew, eh)
    
    def press(self, eb, event):
        pos_list = getattr(self.stream_area, 'pos_gaussians' + self.draw_functions[self.state])
        self.selected_gaussian_index = None
        for i, pos in enumerate(pos_list):
            if ((pos[0] - event.x)**2 + (pos[1] - event.y)**2) <= 20:
                self.selected_gaussian_index = i
        self.drawarea.queue_draw()
            
