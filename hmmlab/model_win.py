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

import os
from os.path import join, abspath, dirname
from gi.repository import Gtk, Gdk, GdkPixbuf
import gtklib

class ModelWindow(gtklib.ObjGetter):
    def __init__(self, model, main_win):
        path = join(dirname(abspath(__file__)), 'glade')
        gtklib.ObjGetter.__init__(self, join(path, "model_win.glade"), self.get_signals())
        self.model = model.model if hasattr(model, 'model') else model
        self.multiplier = 1
        self.image.override_background_color(Gtk.StateFlags.NORMAL, Gdk.RGBA(255,255,255))
        self.main_win = main_win
        self.loaded_state = None
        self.selected_states = []
        self.fill_states_table()
        self.window.show()


    def get_signals(self):
        signals = {"show" : self.show,
                "delete" : self.delete,
                "size_changed" : self.size_changed,
                "scroll" : self.scroll,
                "state_select" : self.state_select,
                "selection_changed" : self.selection_changed,
                "gauss_toggled" : self.gauss_toggled,
                "gauss_select" : self.gauss_select,
                "add_prechod" : self.add_prechod}
        return signals
    
    def delete(self, widget=None, event=None):
        self.main_win.models_windows.remove(self)

    def scroll(self, widget, event):
        self.multiplier = max(0.1, self.multiplier - 0.05*event.delta_y)
        self.load()

    def show(self, widget):
        self.load()
        self.fill_states_table()

    def get_selected_states(self):
        result = []
        for i in range(len(self.model.states)):
            if self.model.modelset.is_selected(self.model, i):
                result.append(i)
        return result

    def load(self):
        alloc = self.image.get_allocation()
        path = self.model.create_image()
        pixbuf = GdkPixbuf.Pixbuf.new_from_file(path)
        width = pixbuf.get_width()
        height = pixbuf.get_height()
        if alloc.width < width:
            height *= alloc.width / width
            width = alloc.width
        if height > alloc.height:
            width = pixbuf.get_width()
            height = pixbuf.get_height()
            width *= alloc.height / height
            height = alloc.height
        pixbuf = GdkPixbuf.Pixbuf.new_from_file_at_scale(path, self.multiplier * width, self.multiplier * height, False)
        self.image.set_from_pixbuf(pixbuf)
        os.remove(path)

    def size_changed(self, widget, alloc):
        self.load()

    def fill_states_table(self, *args):
        self.states_store.clear()
        self.selected_states = self.get_selected_states()
        for i, state in enumerate(self.model.states):
            self.states_store.append([i, state.name, sum([len(s.gaussians) for s in state.streams]), i in self.selected_states])

    def state_select(self, widget, path):
        self.states_store[path][3] = not self.states_store[path][3]
        if self.states_store[path][3]:
            self.model.states[int(path[0])].select_gaussians(True)
        else:
            self.model.states[int(path[0])].unselect_gaussians(True)
        for drawarea in self.main_win.visual_win.streams:
            drawarea.selected_gaussian_index = None
        self.main_win.visual_win.refresh()

    def selection_changed(self, treeview):
        selection = self.states_view.get_selection()
        if selection is not None:
            _, it = selection.get_selected()
            if it is not None:
                state_index = self.states_store.get_value(it, 0)
                self.load_state(self.model.states[state_index])
        else:
            self.name_label.set_text('')
            self.notebook.set_sensitive(False)


    def load_state(self, state):
        for s in self.main_win.visual_win.streams:
            s.select_state(state)
        self.notebook.set_sensitive(True)
        self.loaded_state = state
        self.name_label.set_text('Meno: ' + state.name)
        self.gaussians_store.clear()
        i = 0
        for str_index, stream in enumerate(state.streams):
            for g in stream.gaussians:
                self.gaussians_store.append([i, str_index, g.name, self.model.modelset.is_selected(g)])
                i += 1
        self.prechod_store.clear()
        state_index = self.model.states.index(state)
        for i in range(len(self.model.states)):
            h = self.model.trans_mat(state_index+1, i+1)
            if h > 0:
                self.prechod_store.append([i, self.model.states[i].name, h])

    def gauss_toggled(self, widget, path):
        self.gaussians_store[path][3] = not self.gaussians_store[path][3]
        if self.loaded_state is not None:
            self.loaded_state.get_gaussian(self.gaussians_store[path][0], True)
        self.fill_states_table()
        self.main_win.visual_win.refresh()

    def add_prechod(self, button):
        state_index = self.states_store[self.combobox1.get_active()][0]
        pin = False
        for row in self.prechod_store:
            pin |= row[0] == state_index
        if not pin:
            self.prechod_store.append([state_index, self.model.states[state_index].name, 0.1])
            self.model.trans_mat(self.model.states.index(self.loaded_state)+1, state_index+1, 0.1)
            self.load()

    def gauss_select(self, treeview):
        selection = self.gaussians_view.get_selection()
        if selection is not None:
            _, it = selection.get_selected()
            if it is not None:
                gauss_index = self.gaussians_store.get_value(it, 0)
                str_index = self.gaussians_store.get_value(it, 1)
                if self.loaded_state is not None:
                    gauss = self.loaded_state.get_gaussian(gauss_index, False)
                    if gauss is not None:
                        g_index = -1
                        for i, g in enumerate(self.main_win.visual_win.streams[str_index].stream_area.selected_gaussians):
                            if gauss.name == g.name:
                                g_index = i
                        if g_index != -1:
                            self.main_win.visual_win.streams[str_index].selected_gaussian_index = g_index
                            self.main_win.visual_win.refresh()
