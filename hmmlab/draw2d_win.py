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
from random import random
from os.path import expanduser, join, exists, abspath, dirname
import gtklib

class Draw2DWindow(gtklib.ObjGetter):
    def __init__(self, stream_index, dim1, dim2, modelset, main_win):
        path = join(dirname(abspath(__file__)), 'glade')
        gtklib.ObjGetter.__init__(self, join(path, "draw2d.glade"), self.get_signals())
        self.stream_index = stream_index
        self.dim1 = dim1
        self.dim2 = dim2
        self.modelset = modelset
        self.stream_area = self.modelset.stream_areas[self.stream_index]
        self.data = self.stream_area.get_data_2D(self.dim1, self.dim2)
        self.minx = min([d[0] for d in self.data])
        self.maxx = max([d[0] for d in self.data])
        self.miny = min([d[1] for d in self.data])
        self.maxy = max([d[1] for d in self.data])
        self.window.set_transient_for(main_win)
        self.window.show()

    def get_signals(self):
        signals = {"draw" : self.draw,
                   "scroll" : self.scroll}
        return signals

    def destroy(self, *args):
        self.window.destroy()

    def __del__(self):
        print('Draw2DWindow.__del__')
        self.destroy()

    def scroll(self, *args):
        print(args)

    def draw(self, drawarea, cr):
        #clear
        allocation = drawarea.get_allocation()
        awidth = allocation.width
        aheight = allocation.height
        cr.rectangle(0, 0, allocation.width, allocation.height)
        cr.set_source_rgb(0, 0, 0)
        cr.fill()
        data_width = self.maxx - self.minx
        data_height = self.maxy - self.miny
        xscale = awidth / data_width
        yscale = aheight / data_height
        awidth /= 2
        aheight /= 2

        cr.set_source_rgb (255, 200, 0)
        for d in self.data:
            x = d[0] * xscale + awidth
            y = d[1] * yscale + aheight
            print(x,y)
            cr.move_to(x,y)
            cr.arc(x, y, 3, 0, 2*math.pi);
            cr.fill()
        for gauss in self.stream_area.selected_gaussians:
            mx = gauss.mean[self.dim1] * xscale + awidth
            my = gauss.mean[self.dim2] * yscale + aheight
            r,g,b = random(), random(), random()
            cr.set_source_rgba(r, g, b, 1)
            cr.move_to(mx, my)
            cr.arc(mx, my, 3, 0, 2*math.pi);
            cr.fill()
            vx = math.sqrt(gauss.covariance(self.dim1, self.dim1)) * xscale
            vy = math.sqrt(gauss.covariance(self.dim2, self.dim2)) * yscale
            gtklib.cairo_ellipse(cr, mx - vx / 2, my - vy / 2, vx, vy, r, g, b)
