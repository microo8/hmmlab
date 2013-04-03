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
from gi.repository import Gtk, Gdk
try:
    from hmmlablib import libhmm
    import gtklib
except ImportError:
    from hmmlab.hmmlablib import libhmm
    from hmmlab import gtklib

DATA_DIAMETER = 4
GAUSS_DIAMETER = 5

class DrawArea(gtklib.ObjGetter):
    '''Trieda drawingarea, ktora vykresluje jeden stream vo visualnom okne'''
    def __init__(self, main_win, stream_area):
        '''Vytvori visualne okno'''
        self.main_win = main_win
        self.stream_area = stream_area
        path = join(dirname(abspath(__file__)), 'glade')
        gtklib.ObjGetter.__init__(self, join(path, 'drawarea.glade'), self.get_signals())
        self.selected_train_gaussians = []
        self.selected_train_data = []
        self.selected_gaussian_index = None
        self.selected_state = None
        self.state = 'graphviz'
        self.draw_functions = {'graphviz' : '',
                               'pca' : '_pca',
                               'prob' : '_prob'}

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
        if self.state == 'prierez':
            dim1 = int(self.adjustment1.get_value() - 1)
            dim2 = int(self.adjustment2.get_value() - 1)
            data = self.stream_area.get_data_2D(dim1, dim2)
            if len(data) > 0:
                minx = min([d[0] for d in data])
                maxx = max([d[0] for d in data])
                miny = min([d[1] for d in data])
                maxy = max([d[1] for d in data])
                allocation = drawarea.get_allocation()
                awidth = allocation.width
                aheight = allocation.height
                data_width = maxx - minx
                data_height = maxy - miny
                xscale = float(awidth) / data_width
                yscale = float(aheight) / data_height
                awidth /= 2
                aheight /= 2
                for i, d in enumerate(data):
                    x = d[0] * xscale + awidth - (DATA_DIAMETER / 2.0)
                    y = d[1] * yscale + aheight - (DATA_DIAMETER / 2.0)
                    if self.selected_gaussian_index is not None and self.selected_gaussian_index >= self.stream_area.selected_gaussians.size():
                        self.selected_gaussian_index = None
                    if self.selected_gaussian_index is not None and i in self.stream_area.selected_gaussians[self.selected_gaussian_index].my_data:
                        cr.set_source_rgb(0, 100, 200)
                    elif i in self.selected_train_data:
                        cr.set_source_rgb(255, 255, 255)
                    else:
                        cr.set_source_rgb (255, 200, 0)
                    cr.move_to(x,y)
                    cr.arc(x, y, DATA_DIAMETER, 0, 2*math.pi);
                    cr.fill()
                for i, gauss in enumerate(self.stream_area.selected_gaussians):
                    mx = gauss.mean[dim1] * xscale + awidth - (GAUSS_DIAMETER / 2.0)
                    my = gauss.mean[dim2] * yscale + aheight - (GAUSS_DIAMETER / 2.0)
                    if self.selected_gaussian_index == i:
                        cr.set_source_rgb(0, 0, 255)
                    elif i in self.selected_train_gaussians:
                        cr.set_source_rgb (60, 0, 60)
                    else:
                        cr.set_source_rgb (255, 0, 0)
                    cr.move_to(mx, my)
                    cr.arc(mx, my, GAUSS_DIAMETER, 0, 2*math.pi);
                    cr.fill()
                    vx = math.sqrt(gauss.covariance(dim1, dim1)) * xscale
                    vy = math.sqrt(gauss.covariance(dim2, dim2)) * yscale
                    gtklib.cairo_ellipse(cr, mx - vx / 2, my - vy / 2, vx, vy)
        else:
            pos_list = getattr(self.stream_area, 'pos_data' + self.draw_functions[self.state])
            for i, pos in enumerate(pos_list):
                x = pos[0]
                y = pos[1]
                if self.selected_gaussian_index is not None and self.selected_gaussian_index >= self.stream_area.selected_gaussians.size():
                    self.selected_gaussian_index = None
                if self.selected_gaussian_index is not None and i in self.stream_area.selected_gaussians[self.selected_gaussian_index].my_data:
                    cr.set_source_rgb(0, 100, 200)
                elif i in self.selected_train_data:
                    cr.set_source_rgb(255, 255, 255)
                else:
                    cr.set_source_rgb (255, 200, 0)
                cr.move_to(x,y)
                cr.arc(x, y, 4, 0, 2*math.pi)
                cr.fill()
            pos_list = getattr(self.stream_area, 'pos_gaussians' + self.draw_functions[self.state])
            for i, pos in enumerate(pos_list):
                x = pos[0]
                y = pos[1]
                if self.selected_gaussian_index == i:
                    cr.set_source_rgb(0, 0, 255)
                elif i in self.selected_train_gaussians:
                    cr.set_source_rgb (60, 0, 60)
                else:
                    cr.set_source_rgb (255, 0, 0)
                cr.move_to(x,y)
                cr.arc(x, y, 5, 0, 2*math.pi);
                cr.fill()
                if self.state == 'pca':
                    ew = self.stream_area.pos_gaussians_var_pca[i][0]
                    eh = self.stream_area.pos_gaussians_var_pca[i][1]
                    gtklib.cairo_ellipse(cr, x - ew / 2., y - eh / 2., ew, eh)
    
    def press(self, eb, event):
        if not event.get_state() & Gdk.ModifierType.CONTROL_MASK and event.button == 1:
            self.selected_train_gaussians = []
            self.selected_train_data = []
        if len(self.selected_train_data) > 0 and len(self.selected_train_gaussians) > 0 and event.button == 3:
            self.cluster_popup(event)
        if self.state == 'prierez':
            dim1 = int(self.adjustment1.get_value() - 1)
            dim2 = int(self.adjustment2.get_value() - 1)
            data = self.stream_area.get_data_2D(dim1, dim2)
            if len(data) > 0:
                minx = min([d[0] for d in data])
                maxx = max([d[0] for d in data])
                miny = min([d[1] for d in data])
                maxy = max([d[1] for d in data])
                allocation = self.drawarea.get_allocation()
                awidth = allocation.width
                aheight = allocation.height
                data_width = maxx - minx
                data_height = maxy - miny
                xscale = awidth / data_width
                yscale = aheight / data_height
                awidth /= 2
                aheight /= 2
                for i, gauss in enumerate(self.stream_area.selected_gaussians):
                    mx = gauss.mean[dim1] * xscale + awidth
                    my = gauss.mean[dim2] * yscale + aheight
                    if ((mx - event.x)**2 + (my - event.y)**2) <= GAUSS_DIAMETER**2:
                        if event.get_state() & Gdk.ModifierType.CONTROL_MASK:
                            self.selected_train_gaussians.append(i)
                        else:
                            if event.type == Gdk.EventType.BUTTON_PRESS and event.button == 3:
                                self.open_menu(event, i)
                            elif event.type == Gdk.EventType._2BUTTON_PRESS:
                                self.main_win.open_window_from_gauss(i)
                            else:
                                self.selected_gaussian_index = i
                                self.selected_state = None
                            break
                for i, d in enumerate(self.stream_area.data):
                    dx = d[dim1] * xscale + awidth
                    dy = d[dim2] * yscale + aheight
                    if ((dx - event.x)**2 + (dy - event.y)**2) <= DATA_DIAMETER**2:
                        if event.get_state() & Gdk.ModifierType.CONTROL_MASK:
                            self.selected_train_data.append(i)
                            break
                        else:
                            for j, gauss in enumerate(self.stream_area.selected_gaussians):
                                if i in gauss.my_data:
                                    self.selected_gaussian_index = j
                                    self.selected_state = None
                                    break
                            break
                self.drawarea.queue_draw()
        else:
            pos_list = getattr(self.stream_area, 'pos_gaussians' + self.draw_functions[self.state])
            self.selected_gaussian_index = None
            self.selected_state = None
            for i, pos in enumerate(pos_list):
                if ((pos[0] - event.x)**2 + (pos[1] - event.y)**2) <= GAUSS_DIAMETER**2:
                    if event.get_state() & Gdk.ModifierType.CONTROL_MASK:
                        self.selected_train_gaussians.append(i)
                    else:
                        if event.type == Gdk.EventType.BUTTON_PRESS and event.button == 3:
                            self.open_menu(event, i)
                        elif event.type == Gdk.EventType._2BUTTON_PRESS:
                            gauss = self.stream_area.selected_gaussians[i]
                            self.main_win.open_window_from_gauss(gauss)
                        else:
                            self.selected_gaussian_index = i
                            self.selected_state = None
            pos_list = getattr(self.stream_area, 'pos_data' + self.draw_functions[self.state])
            for i, pos in enumerate(pos_list):
                if ((pos[0] - event.x)**2 + (pos[1] - event.y)**2) <= DATA_DIAMETER**2:
                    if event.get_state() & Gdk.ModifierType.CONTROL_MASK:
                        self.selected_train_data.append(i)
                        break
                    else:
                        for j, gauss in enumerate(self.stream_area.selected_gaussians):
                            if i in gauss.my_data:
                                self.selected_gaussian_index = j
                                self.selected_state = None
                                break
                        break
            self.drawarea.queue_draw()

    def open_menu(self, event, gauss_index):
        menu = Gtk.Menu()
        menu_item = Gtk.MenuItem('Rozdeliť gaussián')
        menu.append(menu_item)
        menu_item.show()
        menu_item.connect("activate", self.menuitem_response, gauss_index)
        menu_item = Gtk.MenuItem('Odstrániť gaussián')
        menu.append(menu_item)
        menu_item.show()
        menu_item.connect("activate", self.delete_gaussian, gauss_index)
        menu.popup(None,
                   None, 
                   lambda menu, data: (event.get_root_coords()[0],
                                       event.get_root_coords()[1],
                                       True),
                   None,
                   event.button,
                   event.time)
    
    def menuitem_response(self, widget, gauss_index):
        gauss = list(self.stream_area.selected_gaussians)[gauss_index]
        gauss.divide()
        self.drawarea.queue_draw()

    def delete_gaussian(self, widget, gauss_index):
        gauss = list(self.stream_area.selected_gaussians)[gauss_index]
        self.stream_area.selected_gaussians.erase(gauss)
        for model in self.stream_area.modelset.get_models_with_gaussian(gauss):
            for s in model.states:
                for st in s.streams:
                    st.gaussians.remove_value(gauss)
        self.stream_area.calc_data_gauss()
        self.stream_area.modelset.reset_pos_gauss()
        self.drawarea.queue_draw()
    
    def cluster_popup(self, event):
        menu = Gtk.Menu()
        menu_item = Gtk.MenuItem('Klastrovať gaussiány')
        menu.append(menu_item)
        menu_item.show()
        menu_item.connect("activate", self.menuitem_cluster)
        menu.popup(None,
                   None, 
                   lambda menu, data: (event.get_root_coords()[0],
                                       event.get_root_coords()[1],
                                       True),
                   None,
                   event.button,
                   event.time)
    
    def menuitem_cluster(self, widget):
        train_gauss = libhmm.ListGaussian()
        for i in self.selected_train_gaussians:
            train_gauss.append(self.stream_area.selected_gaussians[i])
        train_data = libhmm.ListVector()
        for i in self.selected_train_data:
            train_data.append(self.stream_area.get_data(i))
        self.stream_area.modelset.gauss_cluster(train_gauss, train_data)
        self.drawarea.queue_draw()

    def select_state(self, state):
        self.selected_state = state
        self.selected_train_data = [i for i in [self.stream_area.data.index(v) for v in self.stream_area.data if state.in_viterbi_data(v)] if i >= 0]
        self.drawarea.queue_draw()
