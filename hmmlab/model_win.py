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
from threading import Thread
from gi.repository import Gtk, Gdk, GdkPixbuf
try:
    from hmmlablib import libhmm
    import gtklib
except ImportError:
    from hmmlab.hmmlablib import libhmm
    from hmmlab import gtklib
  
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
        self.window.set_title("Model " + self.model.name)
        self.entry1.set_text(self.model.name)
        self.button1.set_sensitive(not self.model.is_joined())
        self.button2.set_sensitive(not self.model.is_joined())
        self.button3.set_sensitive(not self.model.is_joined())
        self.button4.set_sensitive(not self.model.is_joined())
        self.combobox1.set_sensitive(not self.model.is_joined())

        self.states_view.enable_model_drag_source(Gdk.ModifierType.BUTTON1_MASK, [], Gdk.DragAction.COPY)
        self.states_view.drag_source_add_text_targets()
        self.states_view.drag_dest_set(Gtk.DestDefaults.ALL, [], Gdk.DragAction.COPY)
        self.states_view.drag_dest_add_text_targets()
        
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
                "add_prechod" : self.add_prechod,
                "remove_prechod" : self.remove_prechod,
                "name_changed" : self.name_changed,
                "add_state" : self.add_state,
                "remove_state" : self.remove_state,
                "states_data_get" : self.states_data_get,
                "states_data_received" : self.states_data_received,
                "train" : self.train}
        return signals

    def name_changed(self, widget):
        if self.model.modelset.objects_dict.find(self.entry1.get_text()) == self.model.modelset.objects_dict.end():
            self.model.modelset.objects_dict.erase(self.model.name)
            self.model.modelset.objects_dict[self.entry1.get_text()] = self.model
            self.model.name = self.entry1.get_text()
            self.window.set_title("Model " + self.model.name)
            self.main_win.fill_models()
            self.main_win.drawarea.queue_draw()
    
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
        self.combo_store.clear()
        self.selected_states = self.get_selected_states()
        for i, state in enumerate(self.model.states):
            self.states_store.append([i, state.name, sum([len(s.gaussians) for s in state.streams]), i in self.selected_states])
            self.combo_store.append([i, state.name])
        self.combo_store.append([len(self.model.states), "last_%s_state" % self.model.name])

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
                if state_index < self.model.states.size():
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
        for i in range(len(self.model.states)+1):
            h = self.model.trans_mat(state_index+1, i+1)
            if h > 0:
                self.prechod_store.append([i,
                                           self.model.states[i].name
                                           if i < len(self.model.states)
                                           else "last_%s_state" % (self.model.name),
                                           h])

    def gauss_toggled(self, widget, path):
        self.gaussians_store[path][3] = not self.gaussians_store[path][3]
        if self.loaded_state is not None:
            self.loaded_state.get_gaussian(self.gaussians_store[path][0], True)
        self.fill_states_table()
        self.main_win.visual_win.refresh()

    def add_prechod(self, button):
        state_index = self.combo_store[self.combobox1.get_active()][0]
        pin = False
        for row in self.prechod_store:
            pin |= row[0] == state_index
        if not pin:
            self.prechod_store.append([state_index, self.combo_store[self.combobox1.get_active()][1], 0.1])
            self.model.trans_mat(self.model.states.index(self.loaded_state)+1, state_index+1, 0.1)
            self.load()

    def remove_prechod(self, button):
        from_model, from_it = self.states_view.get_selection().get_selected()
        to_model, to_it = self.treeview2.get_selection().get_selected()
        from_index = from_model[from_it][0]
        to_index = to_model[to_it][0]
        self.model.trans_mat(from_index + 1, to_index + 1, 0.0)
        to_model.remove(to_it)
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

    def add_state(self, button):
        s = libhmm.State("%s_state_%d" % (self.model.name, self.model.states.size()+2), self.model.modelset)
        for i, stream in enumerate(s.streams):
            g = libhmm.Gaussian(stream.name + "_gaussian0", self.model.modelset, i, 1)
            dimension = self.model.modelset.streams_distribution[i]
            g.mean = libhmm.SVector(g.name+"_mean", self.model.modelset, dimension, 0.0)
            g.mean.randomize()
            g.covariance = libhmm.SMatrix(g.name+"_covariance", self.model.modelset, dimension, dimension, 0.0)
            g.inv_covariance = libhmm.SMatrix(g.covariance.name+"_inv", self.model.modelset, dimension, dimension, 0.0);
            g.covariance.randomize()
            for i in range(dimension):
                for j in range(dimension):
                    if i != j:
                        g.covariance(i, j, 0.0)
                    else:
                        g.inv_covariance(i, j, 1.0 / g.covariance(i,j))
            g.calc_gconst()
            stream.add_gaussian(g, 1.0)
        self.model.add_state(s)
        self.fill_states_table()
        self.load()

    def remove_state(self, button):
        selection = self.states_view.get_selection()
        if selection is not None:
            _, it = selection.get_selected()
            if it is not None:
                state_index = self.states_store.get_value(it, 0)
                self.model.states[state_index].unselect_gaussians()
                self.model.remove_state(self.model.states[state_index])
                self.fill_states_table()
                self.load()

    def states_data_get(self, widget, drag_context, data, info, time):
        _, it = self.states_view.get_selection().get_selected()
        item = self.states_store.get_value(it, 1)
        data.set_text(str(item), -1)

    def states_data_received(self, widget, drag_context, x, y, data, info, time):
        state_name = data.get_text()
        if state_name not in [s.name for s in self.model.states]:
            state = self.model.modelset.get_state(state_name)
            if state is not None:
                self.model.add_state(state)
                self.fill_states_table()
                self.load()

    def train(self, button):
        iterations = int(self.adjustment2.get_value())
        self.box.set_sensitive(False)
        self.progressbar.set_fraction(0.0)
        self.progressbar.set_text("trénujem")
        self.t = Thread(target=self._train, args=(iterations,))
        self.t.start()

    def _train(self, iterations):
        for i in range(iterations):
            self.model.modelset.train_model(self.model)
            Gdk.threads_enter()
            self.progressbar.set_fraction(float(i+1)/float(iterations))
            self.model.modelset.reset_pos_gauss()
            self.main_win.visual_win.refresh()
            self.main_win.data_table.load_data()
            Gdk.threads_leave()
        Gdk.threads_enter()
        self.progressbar.set_text("")
        self.progressbar.set_fraction(0.0)
        self.box.set_sensitive(True)
        self.fill_states_table()
        self.load()
        self.main_win.visual_win.refresh()
        Gdk.threads_leave()
        del self.t
