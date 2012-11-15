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

from gi.repository import Gtk
import math

class ObjGetter:
    def __init__(self, path, signals):
        self.builder = Gtk.Builder()
        self.builder.add_from_file(path)
        self.builder.connect_signals(signals)

    def __getattr__(self, name):
        return self.builder.get_object(name)


def cairo_rounded_rectangle(cr, x, y, width, height, aspect, corner_radius):
    radius = corner_radius / aspect
    degrees = math.pi / 180.0

    cr.new_sub_path()
    cr.arc(x + width - radius, y + radius, radius, -90 * degrees, 0 * degrees)
    cr.arc(x + width - radius, y + height - radius, radius, 0 * degrees, 90 * degrees)
    cr.arc(x + radius, y + height - radius, radius, 90 * degrees, 180 * degrees)
    cr.arc(x + radius, y + radius, radius, 180 * degrees, 270 * degrees)
    cr.close_path() 

def cairo_ellipse(cr, x, y, width, height, r=-1, g=-1, b=-1):
    if r != -1 and g != -1 and b != -1:
        cr.set_source_rgb(r, g, b)
    cr.save()
    cr.translate(x + width / 2., y + height / 2.)
    cr.scale(1. * (width / 2.), 1. * (height / 2.))
    if r != -1 and g != -1 and b != -1:
        cr.set_source_rgb(r, g, b)
    cr.arc(0., 0., 1., 0., 2 * math.pi)
    cr.restore()
    if r != -1 and g != -1 and b != -1:
        cr.set_source_rgba(r, g, b, 1)
    cr.set_line_width(2)
    cr.stroke()
